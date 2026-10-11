"""Phase 7: does FoldX ΔΔG improve pathogenicity prediction?

Every experiment reuses the Phase 6 data, folds, models, preprocessing and
seed; only the feature set changes. ΔΔG is tested in four paired
comparisons, each "reference" vs "reference + ddg":

1. ``baseline``: conventional + substitution features (the Phase 6 model).
2. ``baseline_structure``: baseline + AlphaFold structure features, to test
   whether ΔΔG adds anything beyond burial and confidence.
3. ``intrinsic``: substitution + structure, without the VEP predictors and
   gnomAD that ClinVar curators use as evidence (circularity-free view).
4. ``substitution``: physicochemistry only, the weakest reference.

For each pair it reports the per-fold differences in every metric, and a
95% confidence interval for the pooled out-of-fold ROC-AUC difference. The
interval comes from a cluster bootstrap that resamples genes (1,000
replicates), because variants within a gene are not independent.

Run with ``python -m biodb.phase7_stability``.
"""

from __future__ import annotations

import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from biodb import config  # noqa: E402
from biodb.common import write_json  # noqa: E402
from biodb.modeling import MODEL_PARAMS, feature_set, roc_points, run_experiment, write_metrics_table  # noqa: E402

OUTPUT_DIR = config.METRICS_DIR / "phase7"
FIGURE_DIR = config.FIGURES_DIR / "phase7"
MODELS_DIR = config.MODELS_DIR / "phase7"
PRIMARY_MODEL = "random_forest"  # best baseline model in Phase 6, fixed before adding ΔΔG
BOOTSTRAP_REPLICATES = 1000

REFERENCES: dict[str, tuple[str, ...]] = {
    "baseline": ("conventional", "substitution"),
    "baseline_structure": ("conventional", "substitution", "structure"),
    "intrinsic": ("substitution", "structure"),
    "substitution": ("substitution",),
}
METRIC_NAMES = ["accuracy", "precision", "recall", "f1", "roc_auc", "pr_auc", "mcc", "brier"]


def gene_bootstrap_auc_difference(oof_ref: pd.DataFrame, oof_new: pd.DataFrame, model: str, seed: int = config.SEED) -> dict[str, float]:
    """Pooled OOF ROC-AUC difference (new - reference) with a gene-level bootstrap 95% CI."""
    y = oof_ref["label"].to_numpy()
    p_ref, p_new = oof_ref[model].to_numpy(), oof_new[model].to_numpy()
    genes = oof_ref["gene"].to_numpy()
    unique = np.unique(genes)
    rows_by_gene = {g: np.flatnonzero(genes == g) for g in unique}
    rng = np.random.default_rng(seed)
    diffs = []
    for _ in range(BOOTSTRAP_REPLICATES):
        sample = np.concatenate([rows_by_gene[g] for g in rng.choice(unique, size=len(unique), replace=True)])
        diffs.append(roc_auc_score(y[sample], p_new[sample]) - roc_auc_score(y[sample], p_ref[sample]))
    diffs = np.array(diffs)
    return {
        "delta_roc_auc": float(roc_auc_score(y, p_new) - roc_auc_score(y, p_ref)),
        "ci95_low": float(np.percentile(diffs, 2.5)),
        "ci95_high": float(np.percentile(diffs, 97.5)),
        "fraction_bootstrap_improved": float((diffs > 0).mean()),
    }


def compare(reference: dict[str, Any], enhanced: dict[str, Any], oof_ref: pd.DataFrame, oof_new: pd.DataFrame) -> dict[str, Any]:
    out = {}
    for model in MODEL_PARAMS:
        ref_folds = reference["results"][model]["per_fold"]
        new_folds = enhanced["results"][model]["per_fold"]
        per_metric = {}
        for name in METRIC_NAMES:
            deltas = [n[name] - r[name] for r, n in zip(ref_folds, new_folds)]
            per_metric[name] = {
                "reference": reference["results"][model]["summary"][name]["mean"],
                "with_ddg": enhanced["results"][model]["summary"][name]["mean"],
                "delta_mean": float(np.mean(deltas)),
                "delta_std": float(np.std(deltas, ddof=1)),
                "folds_improved": int(sum(d < 0 if name == "brier" else d > 0 for d in deltas)),
            }
        importance = enhanced["results"][model]["feature_importance"]
        ranked = list(importance)
        out[model] = {
            "metrics": per_metric,
            "roc_auc_gene_bootstrap": gene_bootstrap_auc_difference(oof_ref, oof_new, model),
            "ddg_importance": importance["ddg"],
            "ddg_importance_rank": ranked.index("ddg") + 1,
            "n_features": len(ranked),
        }
    return out


def plot(comparisons: dict[str, Any], oofs: dict[str, pd.DataFrame]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, len(REFERENCES), figsize=(4 * len(REFERENCES), 4), sharey=True)
    for ax, name in zip(axes, REFERENCES):
        for experiment, style in ((name, "-"), (f"{name}_ddg", "--")):
            oof = oofs[experiment]
            fpr, tpr = roc_points(oof["label"].to_numpy(), oof[PRIMARY_MODEL].to_numpy())
            auc = roc_auc_score(oof["label"], oof[PRIMARY_MODEL])
            ax.plot(fpr, tpr, style, label=f"{experiment} ({auc:.3f})")
        ax.plot([0, 1], [0, 1], color="grey", linestyle=":", linewidth=1)
        ax.set(title=name, xlabel="False positive rate")
        ax.legend(loc="lower right", fontsize=7)
    axes[0].set_ylabel("True positive rate")
    fig.suptitle(f"Out-of-fold ROC, {PRIMARY_MODEL}, with and without ΔΔG")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "roc_with_without_ddg.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4))
    labels, centers, lows, highs = [], [], [], []
    for name in REFERENCES:
        for model in MODEL_PARAMS:
            b = comparisons[name][model]["roc_auc_gene_bootstrap"]
            labels.append(f"{name}\n{model.replace('_', ' ')}")
            centers.append(b["delta_roc_auc"])
            lows.append(b["delta_roc_auc"] - b["ci95_low"])
            highs.append(b["ci95_high"] - b["delta_roc_auc"])
    y = np.arange(len(labels))
    ax.errorbar(centers, y, xerr=[lows, highs], fmt="o", capsize=3)
    ax.axvline(0, color="grey", linestyle=":")
    ax.set_yticks(y, labels, fontsize=7)
    ax.invert_yaxis()
    ax.set(xlabel="Δ ROC-AUC from adding ΔΔG (gene bootstrap 95% CI)")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "delta_auc_ci.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    reports: dict[str, dict[str, Any]] = {}
    for name, groups in REFERENCES.items():
        for experiment, extra in ((name, ()), (f"{name}_ddg", ("stability",))):
            features = feature_set(*groups, *extra)
            reports[experiment] = run_experiment(experiment, features, OUTPUT_DIR, MODELS_DIR)
    oofs = {e: pd.read_csv(OUTPUT_DIR / f"{e}_oof_predictions.csv") for e in reports}

    comparisons = {name: compare(reports[name], reports[f"{name}_ddg"], oofs[name], oofs[f"{name}_ddg"]) for name in REFERENCES}
    write_json(OUTPUT_DIR / "ddg_comparison.json", {
        "primary_model": PRIMARY_MODEL,
        "bootstrap": f"{BOOTSTRAP_REPLICATES} gene-level resamples of pooled out-of-fold predictions",
        "comparisons": comparisons,
    })
    write_metrics_table(list(reports.values()), OUTPUT_DIR / "all_experiments_summary.csv")
    plot(comparisons, oofs)

    print(f"\nAdding ddG ({PRIMARY_MODEL}):")
    for name in REFERENCES:
        c = comparisons[name][PRIMARY_MODEL]
        b = c["roc_auc_gene_bootstrap"]
        print(f"  {name:<20} ROC-AUC {c['metrics']['roc_auc']['reference']:.3f} -> {c['metrics']['roc_auc']['with_ddg']:.3f}  "
              f"delta {b['delta_roc_auc']:+.4f} [{b['ci95_low']:+.4f}, {b['ci95_high']:+.4f}]  "
              f"folds improved {c['metrics']['roc_auc']['folds_improved']}/5  ddg rank {c['ddg_importance_rank']}/{c['n_features']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
