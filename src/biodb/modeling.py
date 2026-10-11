"""Shared model training and evaluation for Phases 6-8.

Every experiment uses the same:

- data: ``data/processed/phase5/ml_dataset.csv``;
- folds: 5-fold ``StratifiedGroupKFold`` grouped by gene (no gene appears in
  both training and test data), created once and saved to
  ``results/metrics/folds.csv``;
- models: logistic regression, random forest and XGBoost with fixed,
  documented hyperparameters (no tuning, so no tuning leakage);
- preprocessing: median imputation (and scaling for logistic regression)
  fitted inside each training fold. XGBoost handles missing values natively.

Only the feature set changes between experiments.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score, average_precision_score, brier_score_loss, confusion_matrix, f1_score,
    matthews_corrcoef, precision_score, recall_score, roc_auc_score, roc_curve,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from xgboost import XGBClassifier

from biodb import config
from biodb.common import write_json

N_FOLDS = 5
THRESHOLD = 0.5
DATASET = config.PHASE5_DIR / "ml_dataset.csv"
METADATA = config.PHASE5_DIR / "feature_metadata.json"
FOLDS = config.METRICS_DIR / "folds.csv"

MODEL_PARAMS: dict[str, dict[str, Any]] = {
    "logistic_regression": {"C": 1.0, "max_iter": 5000},
    "random_forest": {"n_estimators": 500, "min_samples_leaf": 5, "max_features": "sqrt", "n_jobs": 4},
    "xgboost": {
        "n_estimators": 400, "max_depth": 4, "learning_rate": 0.05, "subsample": 0.8,
        "colsample_bytree": 0.8, "min_child_weight": 5, "eval_metric": "logloss", "n_jobs": 4,
    },
}


def feature_groups() -> dict[str, list[str]]:
    return json.loads(METADATA.read_text(encoding="utf-8"))["feature_groups"]


def feature_set(*groups: str) -> list[str]:
    available = feature_groups()
    return [f for g in groups for f in available[g]]


def load_dataset() -> pd.DataFrame:
    return pd.read_csv(DATASET)


def make_model(name: str) -> Any:
    params = MODEL_PARAMS[name]
    if name == "logistic_regression":
        return Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", LogisticRegression(random_state=config.SEED, **params)),
        ])
    if name == "random_forest":
        return Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("model", RandomForestClassifier(random_state=config.SEED, **params)),
        ])
    if name == "xgboost":
        return XGBClassifier(random_state=config.SEED, **params)
    raise ValueError(name)


def get_folds(data: pd.DataFrame) -> pd.Series:
    """Fold index per variant. Created once, then always read from disk."""
    if FOLDS.exists():
        saved = pd.read_csv(FOLDS).set_index("variant_id")["fold"]
        if set(saved.index) != set(data["variant_id"]):
            raise ValueError("saved folds do not match the dataset; delete folds.csv only if the dataset changed on purpose")
        return data["variant_id"].map(saved)
    splitter = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=config.SEED)
    folds = pd.Series(-1, index=data.index)
    for fold, (_, test) in enumerate(splitter.split(data, data["label"], groups=data["gene"])):
        folds.iloc[test] = fold
    FOLDS.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"variant_id": data["variant_id"], "gene": data["gene"], "fold": folds}).to_csv(FOLDS, index=False, lineterminator="\n")
    return folds


def metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    predicted = (probability >= THRESHOLD).astype(int)
    return {
        "accuracy": accuracy_score(y, predicted),
        "precision": precision_score(y, predicted, zero_division=0),
        "recall": recall_score(y, predicted),
        "f1": f1_score(y, predicted),
        "roc_auc": roc_auc_score(y, probability),
        "pr_auc": average_precision_score(y, probability),
        "mcc": matthews_corrcoef(y, predicted),
        "brier": brier_score_loss(y, probability),
    }


def cross_validate(data: pd.DataFrame, features: list[str], model_name: str, folds: pd.Series) -> dict[str, Any]:
    """Out-of-fold predictions and per-fold metrics for one model and feature set."""
    X, y = data[features], data["label"].to_numpy()
    oof = np.zeros(len(data))
    per_fold = []
    for fold in range(N_FOLDS):
        train, test = (folds != fold).to_numpy(), (folds == fold).to_numpy()
        assert not set(data.loc[train, "gene"]) & set(data.loc[test, "gene"]), "gene leaked across folds"
        model = make_model(model_name)
        model.fit(X[train], y[train])
        oof[test] = model.predict_proba(X[test])[:, 1]
        per_fold.append({"fold": fold, "n_test": int(test.sum()), **metrics(y[test], oof[test])})
    summary = {
        name: {"mean": float(np.mean([f[name] for f in per_fold])), "std": float(np.std([f[name] for f in per_fold], ddof=1))}
        for name in per_fold[0] if name not in ("fold", "n_test")
    }
    tn, fp, fn, tp = confusion_matrix(y, (oof >= THRESHOLD).astype(int)).ravel()
    return {
        "model": model_name,
        "features": features,
        "per_fold": per_fold,
        "summary": summary,
        "pooled_out_of_fold": metrics(y, oof),
        "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
        "oof": oof,
    }


def feature_importance(model: Any, features: list[str]) -> dict[str, float]:
    estimator = model.named_steps["model"] if isinstance(model, Pipeline) else model
    if hasattr(estimator, "feature_importances_"):
        values = estimator.feature_importances_
    else:
        values = np.abs(estimator.coef_[0])
    return dict(sorted(zip(features, map(float, values)), key=lambda kv: -kv[1]))


def run_experiment(
    name: str,
    features: list[str],
    output_dir: Path,
    models_dir: Path,
) -> dict[str, Any]:
    """Cross-validate all models on one feature set, fit final models on all data, save everything."""
    data = load_dataset()
    folds = get_folds(data)
    output_dir.mkdir(parents=True, exist_ok=True)
    models_dir.mkdir(parents=True, exist_ok=True)

    results = {}
    oof_columns = {"variant_id": data["variant_id"], "gene": data["gene"], "label": data["label"], "fold": folds}
    for model_name in MODEL_PARAMS:
        result = cross_validate(data, features, model_name, folds)
        oof_columns[model_name] = result.pop("oof")
        final = make_model(model_name).fit(data[features], data["label"])
        joblib.dump({"model": final, "features": features, "experiment": name}, models_dir / f"{name}_{model_name}.joblib")
        result["feature_importance"] = feature_importance(final, features)
        results[model_name] = result
        s = result["summary"]
        print(f"  {name:<28} {model_name:<20} ROC-AUC {s['roc_auc']['mean']:.3f}±{s['roc_auc']['std']:.3f}  "
              f"F1 {s['f1']['mean']:.3f}  MCC {s['mcc']['mean']:.3f}", flush=True)

    pd.DataFrame(oof_columns).to_csv(output_dir / f"{name}_oof_predictions.csv", index=False, lineterminator="\n")
    best = max(results, key=lambda m: results[m]["summary"]["roc_auc"]["mean"])
    report = {
        "experiment": name,
        "features": features,
        "n_features": len(features),
        "cv": f"{N_FOLDS}-fold StratifiedGroupKFold by gene, seed {config.SEED}",
        "threshold": THRESHOLD,
        "model_params": MODEL_PARAMS,
        "best_model_by_roc_auc": best,
        "results": results,
    }
    write_json(output_dir / f"{name}_metrics.json", report)
    return report


def roc_points(y: np.ndarray, probability: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fpr, tpr, _ = roc_curve(y, probability)
    return fpr, tpr


def write_metrics_table(reports: list[dict[str, Any]], path: Path) -> None:
    """CSV of mean±std CV metrics for every experiment and model."""
    names = ["accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc", "mcc", "brier"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["experiment", "model", "n_features"] + [f"{n}_{s}" for n in names for s in ("mean", "std")])
        for report in reports:
            for model_name, result in report["results"].items():
                s = result["summary"]
                writer.writerow([report["experiment"], model_name, report["n_features"]] +
                                [f"{s[n][k]:.4f}" for n in names for k in ("mean", "std")])
