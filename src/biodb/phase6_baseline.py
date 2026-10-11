"""Phase 6: baseline pathogenicity model on conventional features only.

The baseline uses the ``conventional`` (VEP scores, conservation, gnomAD) and
``substitution`` (amino-acid physicochemistry) feature groups. It deliberately
excludes the AlphaFold structure features and FoldX ΔΔG, which Phase 7 adds.

Run with ``python -m biodb.phase6_baseline``.
"""

from __future__ import annotations

import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from biodb import config  # noqa: E402
from biodb.modeling import feature_set, roc_points, run_experiment, write_metrics_table  # noqa: E402

EXPERIMENT = "baseline"
OUTPUT_DIR = config.METRICS_DIR / "phase6"
FIGURE_DIR = config.FIGURES_DIR / "phase6"
MODELS_DIR = config.MODELS_DIR / "phase6"


def plot(report: dict, oof: pd.DataFrame) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    y = oof["label"].to_numpy()

    fig, ax = plt.subplots(figsize=(5, 5))
    for model_name, result in report["results"].items():
        fpr, tpr = roc_points(y, oof[model_name].to_numpy())
        ax.plot(fpr, tpr, label=f"{model_name} (AUC {result['summary']['roc_auc']['mean']:.3f})")
    ax.plot([0, 1], [0, 1], color="grey", linestyle=":", linewidth=1)
    ax.set(xlabel="False positive rate", ylabel="True positive rate", title="Baseline, out-of-fold ROC (gene-grouped CV)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "roc_baseline.png", dpi=150)
    plt.close(fig)

    best = report["best_model_by_roc_auc"]
    cm = report["results"][best]["confusion_matrix"]
    matrix = np.array([[cm["tn"], cm["fp"]], [cm["fn"], cm["tp"]]])
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.imshow(matrix, cmap="Blues")
    for (i, j), value in np.ndenumerate(matrix):
        ax.text(j, i, str(value), ha="center", va="center", color="white" if value > matrix.max() / 2 else "black")
    ax.set(xticks=[0, 1], yticks=[0, 1], xticklabels=["benign", "pathogenic"], yticklabels=["benign", "pathogenic"],
           xlabel="Predicted", ylabel="True", title=f"Baseline {best}, threshold 0.5")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "confusion_baseline.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    features = feature_set("conventional", "substitution")
    print(f"Baseline: {len(features)} features")
    report = run_experiment(EXPERIMENT, features, OUTPUT_DIR, MODELS_DIR)
    write_metrics_table([report], OUTPUT_DIR / "baseline_summary.csv")
    plot(report, pd.read_csv(OUTPUT_DIR / f"{EXPERIMENT}_oof_predictions.csv"))
    print(f"Best by ROC-AUC: {report['best_model_by_roc_auc']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
