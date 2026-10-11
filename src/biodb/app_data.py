"""Data access for the Streamlit app (Phase 9).

The app never trains anything. It reads versioned project outputs:

- ``data/processed/phase5/ml_dataset.csv``: features per variant;
- ``data/processed/phase1/cohort.csv``: ClinVar details per variant;
- ``data/processed/phase1/selected_genes.csv``: UniProt and AlphaFold IDs;
- ``results/metrics/phase7/{baseline,baseline_ddg}_oof_predictions.csv``:
  out-of-fold probabilities from the gene-grouped cross-validation. Each
  probability was produced by a model that never saw the variant's gene.
- ``models/phase7/*.joblib``: the final models refitted on all variants.

Structures are fetched from AlphaFold DB on demand rather than shipped.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Any

import joblib
import pandas as pd
import requests

from biodb import config

MODEL = "xgboost"  # the only tree model small enough to version (see Phase 6)
EXPERIMENTS = {"baseline": "Baseline (VEP + substitution)", "baseline_ddg": "Baseline + FoldX ΔΔG"}
THRESHOLD = 0.5


@dataclass(frozen=True)
class Prediction:
    experiment: str
    label: str
    out_of_fold: float
    final_model: float

    @property
    def predicted_class(self) -> str:
        return "Pathogenic" if self.out_of_fold >= THRESHOLD else "Benign"


@lru_cache(maxsize=1)
def tables() -> dict[str, pd.DataFrame]:
    dataset = pd.read_csv(config.PHASE5_DIR / "ml_dataset.csv")
    cohort = pd.read_csv(config.PHASE1_DIR / "cohort.csv")
    genes = pd.read_csv(config.PHASE1_DIR / "selected_genes.csv")
    oof = {
        name: pd.read_csv(config.METRICS_DIR / "phase7" / f"{name}_oof_predictions.csv").set_index("variant_id")
        for name in EXPERIMENTS
    }
    return {"dataset": dataset.set_index("variant_id"), "cohort": cohort.set_index("variant_id"),
            "genes": genes.set_index("gene"), **{f"oof_{k}": v for k, v in oof.items()}}


@lru_cache(maxsize=None)
def final_model(experiment: str) -> dict[str, Any]:
    return joblib.load(config.MODELS_DIR / "phase7" / f"{experiment}_{MODEL}.joblib")


def gene_list() -> list[str]:
    return sorted(tables()["genes"].index)


def variants_for(gene: str) -> list[str]:
    data = tables()["dataset"]
    subset = data[data["gene"] == gene].sort_values("protein_position")
    return list(subset.index)


def parse_variant(gene: str, text: str) -> str | None:
    """Accept 'R175H' or 'p.R175H' (one-letter codes) and return the cohort variant_id, if present."""
    change = text.strip().removeprefix("p.").upper()
    variant_id = f"{gene}:p.{change}"
    return variant_id if variant_id in tables()["dataset"].index else None


def variant_details(variant_id: str) -> dict[str, Any]:
    t = tables()
    row = t["dataset"].loc[variant_id]
    clinvar = t["cohort"].loc[variant_id]
    gene = t["genes"].loc[row["gene"]]
    return {
        "row": row,
        "clinvar_significance": clinvar["ClinicalSignificance"],
        "review_status": clinvar["ReviewStatus"],
        "variation_id": int(clinvar["VariationID"]),
        "hgvs": clinvar["Name"],
        "uniprot": gene["uniprot_accession"],
        "alphafold_entry": gene["alphafold_entry"],
        "cif_url": gene["alphafold_cif_url"],
    }


def predictions(variant_id: str) -> list[Prediction]:
    t = tables()
    row = t["dataset"].loc[[variant_id]]
    out = []
    for experiment, label in EXPERIMENTS.items():
        bundle = final_model(experiment)
        final = float(bundle["model"].predict_proba(row[bundle["features"]])[0, 1])
        out.append(Prediction(experiment, label, float(t[f"oof_{experiment}"].loc[variant_id, MODEL]), final))
    return out


@lru_cache(maxsize=16)
def fetch_structure(url: str) -> str:
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text
