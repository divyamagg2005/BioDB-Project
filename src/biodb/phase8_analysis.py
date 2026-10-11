"""Phase 8: scientific analysis of the ΔΔG results.

Uses only saved outputs from Phases 5-7, no retraining:

- single-feature ROC-AUC of every feature (how strong each signal is alone);
- ΔΔG discrimination within AlphaFold-confidence and burial strata;
- Spearman correlation of ΔΔG with the structural features (redundancy);
- per-gene ΔΔG ROC-AUC (is the signal consistent across proteins?);
- where ΔΔG could help most: out-of-fold performance with and without ΔΔG on
  variants the baseline is unsure about, and on variants whose AlphaMissense
  class is "ambiguous";
- figures: ΔΔG distributions, feature importance, final comparison.

Run with ``python -m biodb.phase8_analysis``.
"""

from __future__ import annotations

import json
import sys
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from biodb import config  # noqa: E402
from biodb.common import write_json  # noqa: E402
from biodb.modeling import feature_groups, load_dataset  # noqa: E402

PHASE7 = config.METRICS_DIR / "phase7"
OUTPUT_DIR = config.METRICS_DIR / "phase8"
FIGURE_DIR = config.FIGURES_DIR / "phase8"
MODEL = "random_forest"
# AlphaMissense class boundaries (Cheng et al. 2023): < 0.34 benign, > 0.564 pathogenic.
AM_AMBIGUOUS = (0.34, 0.564)


def auc(y: pd.Series, score: pd.Series) -> float | None:
    mask = score.notna()
    if y[mask].nunique() < 2:
        return None
    return float(roc_auc_score(y[mask], score[mask]))


def single_feature_auc(data: pd.DataFrame) -> dict[str, dict[str, Any]]:
    """ROC-AUC of each feature alone, oriented so that > 0.5 means informative."""
    out = {}
    for group, features in feature_groups().items():
        for f in features:
            raw = auc(data["label"], data[f])
            if raw is None:
                continue
            out[f] = {"group": group, "auc": max(raw, 1 - raw), "direction": "higher = pathogenic" if raw >= 0.5 else "lower = pathogenic"}
    return dict(sorted(out.items(), key=lambda kv: -kv[1]["auc"]))


def strata(data: pd.DataFrame) -> dict[str, Any]:
    rows = {
        "pLDDT >= 70": data["plddt"] >= 70,
        "pLDDT < 70": data["plddt"] < 70,
        "buried (RSA < 0.25)": data["rsa"] < 0.25,
        "exposed (RSA >= 0.25)": data["rsa"] >= 0.25,
        "pLDDT >= 70 and buried": (data["plddt"] >= 70) & (data["rsa"] < 0.25),
        "pLDDT >= 70 and exposed": (data["plddt"] >= 70) & (data["rsa"] >= 0.25),
    }
    out = {}
    for name, mask in rows.items():
        sub = data[mask]
        out[name] = {
            "n": int(len(sub)),
            "pathogenic": int(sub["label"].sum()),
            "ddg_auc": auc(sub["label"], sub["ddg"]),
            "median_ddg_pathogenic": float(sub.loc[sub["label"] == 1, "ddg"].median()),
            "median_ddg_benign": float(sub.loc[sub["label"] == 0, "ddg"].median()),
        }
    return out


def redundancy(data: pd.DataFrame) -> dict[str, float]:
    return {
        f: float(spearmanr(data["ddg"], data[f], nan_policy="omit").statistic)
        for f in ["plddt", "plddt_window", "rsa", "neighbor_count", "blosum62", "abs_volume_change",
                  "alphamissense", "revel", "cadd_phred"]
    }


def per_gene(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene, sub in data.groupby("gene"):
        rows.append({
            "gene": gene,
            "ddg_auc": auc(sub["label"], sub["ddg"]),
            "median_plddt": float(sub["plddt"].median()),
            "median_ddg_pathogenic": float(sub.loc[sub["label"] == 1, "ddg"].median()),
            "median_ddg_benign": float(sub.loc[sub["label"] == 0, "ddg"].median()),
        })
    return pd.DataFrame(rows).sort_values("ddg_auc", ascending=False)


def hard_cases(data: pd.DataFrame) -> dict[str, Any]:
    """With vs without ΔΔG on the variants where the baseline is least certain."""
    out = {}
    for reference in ("baseline", "baseline_structure"):
        ref = pd.read_csv(PHASE7 / f"{reference}_oof_predictions.csv").set_index("variant_id")
        new = pd.read_csv(PHASE7 / f"{reference}_ddg_oof_predictions.csv").set_index("variant_id")
        merged = data.set_index("variant_id").join(ref[[MODEL]].rename(columns={MODEL: "p_ref"})).join(
            new[[MODEL]].rename(columns={MODEL: "p_new"}))
        subsets = {
            "all": pd.Series(True, index=merged.index),
            "baseline uncertain (0.2 <= p <= 0.8)": merged["p_ref"].between(0.2, 0.8),
            "AlphaMissense ambiguous": merged["alphamissense"].between(*AM_AMBIGUOUS),
            "AlphaMissense missing": merged["alphamissense"].isna(),
        }
        out[reference] = {}
        for name, mask in subsets.items():
            sub = merged[mask]
            out[reference][name] = {
                "n": int(len(sub)),
                "auc_without_ddg": auc(sub["label"], sub["p_ref"]),
                "auc_with_ddg": auc(sub["label"], sub["p_new"]),
                "accuracy_without_ddg": float(((sub["p_ref"] >= 0.5) == sub["label"]).mean()),
                "accuracy_with_ddg": float(((sub["p_new"] >= 0.5) == sub["label"]).mean()),
            }
    return out


def figures(data: pd.DataFrame, genes: pd.DataFrame, single: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    clipped = data.assign(ddg=data["ddg"].clip(-3, 15))
    for ax, (title, mask) in zip(axes, (("pLDDT >= 70", data["plddt"] >= 70), ("pLDDT < 70", data["plddt"] < 70))):
        sub = clipped[mask]
        parts = [sub.loc[sub["label"] == 0, "ddg"], sub.loc[sub["label"] == 1, "ddg"]]
        ax.boxplot(parts, showfliers=False)
        ax.set_xticks([1, 2], [f"benign (n={len(parts[0])})", f"pathogenic (n={len(parts[1])})"])
        ax.axhline(0, color="grey", linestyle=":")
        ax.set_title(f"{title}, ΔΔG AUC {auc(data[mask]['label'], data[mask]['ddg']):.3f}")
    axes[0].set_ylabel("FoldX ΔΔG (kcal/mol, clipped to [-3, 15])")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "ddg_by_label_and_confidence.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 6))
    names = list(single)[::-1]
    colors = {"conventional": "#4c72b0", "substitution": "#8c8c8c", "structure": "#55a868", "stability": "#c44e52"}
    ax.barh(names, [single[n]["auc"] for n in names], color=[colors[single[n]["group"]] for n in names])
    ax.axvline(0.5, color="grey", linestyle=":")
    ax.set(xlim=(0.45, 1.0), xlabel="ROC-AUC of the feature alone", title="Single-feature discrimination")
    ax.legend(handles=[Patch(color=c, label=g) for g, c in colors.items()], loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "single_feature_auc.png", dpi=150)
    plt.close(fig)

    metrics = json.loads((PHASE7 / "baseline_structure_ddg_metrics.json").read_text(encoding="utf-8"))
    importance = metrics["results"][MODEL]["feature_importance"]
    top = list(importance.items())[:15][::-1]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.barh([k for k, _ in top], [v for _, v in top],
            color=["#c44e52" if k == "ddg" else "#4c72b0" for k, _ in top])
    ax.set(xlabel="Random forest impurity importance", title="Top 15 features, baseline + structure + ΔΔG")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "feature_importance_full_model.png", dpi=150)
    plt.close(fig)

    valid = genes.dropna(subset=["ddg_auc"])
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(valid["median_plddt"], valid["ddg_auc"], s=15)
    ax.axhline(0.5, color="grey", linestyle=":")
    ax.set(xlabel="Median pLDDT of the gene's cohort variants", ylabel="ΔΔG ROC-AUC within gene",
           title="Per-gene ΔΔG discrimination (25 + 25 variants each)")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "per_gene_ddg_auc.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    data = load_dataset()
    single = single_feature_auc(data)
    genes = per_gene(data)
    result = {
        "single_feature_auc": single,
        "ddg_strata": strata(data),
        "ddg_spearman_with": redundancy(data),
        "per_gene_ddg_auc_summary": {
            "genes": int(genes["ddg_auc"].notna().sum()),
            "median": float(genes["ddg_auc"].median()),
            "above_0.7": int((genes["ddg_auc"] > 0.7).sum()),
            "below_0.5": int((genes["ddg_auc"] < 0.5).sum()),
            "spearman_with_median_plddt": float(spearmanr(genes["median_plddt"], genes["ddg_auc"], nan_policy="omit").statistic),
        },
        "hard_cases": hard_cases(data),
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    genes.to_csv(OUTPUT_DIR / "per_gene_ddg_auc.csv", index=False, lineterminator="\n", float_format="%.4f")
    write_json(OUTPUT_DIR / "analysis.json", result)
    figures(data, genes, single)

    print("Single-feature AUC (top 8):", {k: round(v["auc"], 3) for k, v in list(single.items())[:8]})
    print("ddg AUC by stratum:", {k: (v["n"], round(v["ddg_auc"], 3)) for k, v in result["ddg_strata"].items()})
    print("ddg Spearman:", {k: round(v, 2) for k, v in result["ddg_spearman_with"].items()})
    print("per gene:", result["per_gene_ddg_auc_summary"])
    for ref, subsets in result["hard_cases"].items():
        for name, v in subsets.items():
            print(f"  {ref:<20} {name:<40} n={v['n']:<5} AUC {v['auc_without_ddg']:.3f} -> {v['auc_with_ddg']:.3f}  "
                  f"acc {v['accuracy_without_ddg']:.3f} -> {v['accuracy_with_ddg']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
