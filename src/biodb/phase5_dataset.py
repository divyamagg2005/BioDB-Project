"""Phase 5: merge all phases into one modelling table.

One row per cohort variant with the label, the gene (used only to group
cross-validation folds) and four feature groups:

- ``conventional``: VEP scores, conservation and population frequency.
- ``substitution``: physicochemical change of the amino-acid substitution.
- ``structure``: AlphaFold-derived residue environment.
- ``stability``: FoldX ΔΔG.

Biological missing values stay missing (empty). Imputation belongs inside
the cross-validated model pipeline, so it is fitted on training folds only.
The one exception is ``gnomad_absent``: absence from gnomAD is itself an
observation, so it is recorded as an explicit 0/1 column here.

Run with ``python -m biodb.phase5_dataset``.
"""

from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from Bio.Align import substitution_matrices

from biodb import config
from biodb.common import sha256_file, write_json

# Kyte & Doolittle 1982 hydropathy.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
# Zamyatnin 1972 residue volumes (A^3).
VOLUME = {
    "G": 60.1, "A": 88.6, "S": 89.0, "C": 108.5, "D": 111.1, "P": 112.7, "N": 114.1, "T": 116.1, "E": 138.4,
    "V": 140.0, "Q": 143.8, "H": 153.2, "M": 162.9, "I": 166.7, "L": 166.7, "K": 168.6, "R": 173.4, "F": 189.9,
    "Y": 193.6, "W": 227.8,
}
# Side-chain charge at physiological pH (histidine treated as neutral).
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1}
BLOSUM62 = substitution_matrices.load("BLOSUM62")

FEATURE_GROUPS: dict[str, list[str]] = {
    "conventional": [
        "cadd_phred", "revel", "alphamissense", "sift_score", "polyphen_score", "gerp_conservation",
        "gnomad_af_max", "gnomad_absent",
    ],
    "substitution": [
        "blosum62", "hydropathy_change", "volume_change", "abs_volume_change", "charge_change",
        "to_proline", "from_proline", "to_glycine", "from_glycine", "from_cysteine", "to_cysteine",
    ],
    "structure": ["plddt", "plddt_window", "rsa", "neighbor_count", "ss_helix", "ss_strand"],
    "stability": ["ddg"],
}
FEATURE_DESCRIPTIONS = {
    "cadd_phred": "CADD v1.7 PHRED-scaled deleteriousness (VEP)",
    "revel": "REVEL ensemble missense score, 0-1 (VEP)",
    "alphamissense": "AlphaMissense pathogenicity, 0-1 (VEP)",
    "sift_score": "SIFT score, 0-1, low = deleterious (VEP)",
    "polyphen_score": "PolyPhen-2 HumVar score, 0-1 (VEP)",
    "gerp_conservation": "GERP++ conservation score (VEP Conservation plugin)",
    "gnomad_af_max": "Max of gnomAD exome/genome allele frequency for the exact allele; empty if absent",
    "gnomad_absent": "1 if the allele has no gnomAD frequency, else 0",
    "blosum62": "BLOSUM62 substitution score (wild type to mutant)",
    "hydropathy_change": "Kyte-Doolittle hydropathy, mutant minus wild type",
    "volume_change": "Residue volume (A^3, Zamyatnin), mutant minus wild type",
    "abs_volume_change": "Absolute residue volume change",
    "charge_change": "Side-chain charge, mutant minus wild type (D/E -1, K/R +1)",
    "to_proline": "1 if the mutant residue is proline",
    "from_proline": "1 if the wild-type residue is proline",
    "to_glycine": "1 if the mutant residue is glycine",
    "from_glycine": "1 if the wild-type residue is glycine",
    "from_cysteine": "1 if the wild-type residue is cysteine",
    "to_cysteine": "1 if the mutant residue is cysteine",
    "plddt": "AlphaFold per-residue confidence, 0-100",
    "plddt_window": "Mean pLDDT over residues i-5..i+5",
    "rsa": "Relative solvent accessibility (Shrake-Rupley SASA / Tien 2013 max)",
    "neighbor_count": "Residues with CB within 10 A of this residue's CB",
    "ss_helix": "1 if the residue is in a helix (DSSP, AlphaFold mmCIF)",
    "ss_strand": "1 if the residue is in a strand",
    "ddg": "FoldX 5.1 BuildModel ΔΔG, kcal/mol, mean of 3 runs (positive = destabilising)",
}
ID_COLUMNS = ["variant_id", "gene", "uniprot_accession", "protein_position", "wt_aa", "mut_aa", "label"]
ALL_FEATURES = [f for group in FEATURE_GROUPS.values() for f in group]


def substitution_features(wt: str, mut: str) -> dict[str, float]:
    return {
        "blosum62": float(BLOSUM62[wt][mut]),
        "hydropathy_change": round(HYDROPATHY[mut] - HYDROPATHY[wt], 2),
        "volume_change": round(VOLUME[mut] - VOLUME[wt], 1),
        "abs_volume_change": round(abs(VOLUME[mut] - VOLUME[wt]), 1),
        "charge_change": CHARGE.get(mut, 0) - CHARGE.get(wt, 0),
        "to_proline": int(mut == "P"),
        "from_proline": int(wt == "P"),
        "to_glycine": int(mut == "G"),
        "from_glycine": int(wt == "G"),
        "from_cysteine": int(wt == "C"),
        "to_cysteine": int(mut == "C"),
    }


def _float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _index(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    index = {r["variant_id"]: r for r in rows}
    if len(index) != len(rows):
        raise ValueError(f"duplicate variant_id in {path}")
    return index


def build_row(cohort: dict[str, str], vep: dict[str, str], structure: dict[str, str], foldx: dict[str, str]) -> dict[str, Any]:
    wt, mut = cohort["ref_aa"], cohort["alt_aa"]
    af = _float(vep.get("gnomad_af_max"))
    row: dict[str, Any] = {
        "variant_id": cohort["variant_id"],
        "gene": cohort["gene"],
        "uniprot_accession": cohort["uniprot_accession"],
        "protein_position": int(cohort["protein_position"]),
        "wt_aa": wt,
        "mut_aa": mut,
        "label": int(cohort["label"]),
        "cadd_phred": _float(vep.get("cadd_phred")),
        "revel": _float(vep.get("revel")),
        "alphamissense": _float(vep.get("alphamissense")),
        "sift_score": _float(vep.get("sift_score")),
        "polyphen_score": _float(vep.get("polyphen_score")),
        "gerp_conservation": _float(vep.get("gerp_conservation")),
        "gnomad_af_max": af,
        "gnomad_absent": int(af is None),
        **substitution_features(wt, mut),
        "plddt": _float(structure.get("plddt")),
        "plddt_window": _float(structure.get("plddt_window")),
        "rsa": _float(structure.get("rsa")),
        "neighbor_count": _float(structure.get("neighbor_count")),
        "ss_helix": int(structure.get("secondary_structure") == "helix"),
        "ss_strand": int(structure.get("secondary_structure") == "strand"),
        "ddg": _float(foldx.get("ddg")) if foldx.get("status") == "ok" else None,
    }
    return row


def validate(rows: list[dict[str, Any]], cohort: list[dict[str, str]], target_genes: int) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: Any = "") -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    labels = Counter(r["label"] for r in rows)
    check("row_per_cohort_variant", len(rows) == len(cohort), len(rows))
    check("unique_variant_id", len({r["variant_id"] for r in rows}) == len(rows))
    check("gene_count", len({r["gene"] for r in rows}) == target_genes, len({r["gene"] for r in rows}))
    check("class_balance", labels[1] == labels[0], dict(labels))
    check("no_label_derived_features", not {"label", "ClinicalSignificance", "ReviewStatus"} & set(ALL_FEATURES))
    check("features_numeric_or_missing", all(
        r[f] is None or isinstance(r[f], (int, float)) for r in rows for f in ALL_FEATURES
    ))
    check("ddg_present_for_all", all(r["ddg"] is not None for r in rows))
    check("gnomad_absent_consistent", all((r["gnomad_af_max"] is None) == bool(r["gnomad_absent"]) for r in rows))
    return checks


def run(output_dir: Path = config.PHASE5_DIR) -> dict[str, Any]:
    cohort_path = config.PHASE1_DIR / "cohort.csv"
    with cohort_path.open(encoding="utf-8", newline="") as handle:
        cohort = list(csv.DictReader(handle))
    vep = _index(config.PHASE2_DIR / "vep_annotations.csv")
    structure = _index(config.PHASE3_DIR / "variant_structure_features.csv")
    foldx = _index(config.PHASE4_DIR / "foldx_ddg.csv")

    rows = [build_row(c, vep.get(c["variant_id"], {}), structure.get(c["variant_id"], {}), foldx.get(c["variant_id"], {}))
            for c in cohort]
    merge_losses = {
        name: sum(c["variant_id"] not in source for c in cohort)
        for name, source in (("vep", vep), ("structure", structure), ("foldx", foldx))
    }
    checks = validate(rows, cohort, config.TARGET_GENES)
    checks.append({"check": "no_merge_losses", "passed": not any(merge_losses.values()), "detail": merge_losses})
    passed = all(c["passed"] for c in checks)

    output_dir.mkdir(parents=True, exist_ok=True)
    dataset_path = output_dir / "ml_dataset.csv"
    with dataset_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ID_COLUMNS + ALL_FEATURES, lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: ("" if v is None else v) for k, v in r.items()} for r in rows)

    missingness = {
        f: {
            "missing": sum(r[f] is None for r in rows),
            "missing_pathogenic": sum(r[f] is None for r in rows if r["label"] == 1),
            "missing_benign": sum(r[f] is None for r in rows if r["label"] == 0),
        }
        for f in ALL_FEATURES
    }
    write_json(output_dir / "feature_metadata.json", {
        "id_columns": ID_COLUMNS,
        "target": "label (1 = pathogenic/likely pathogenic, 0 = benign/likely benign)",
        "group_column": "gene (for gene-grouped cross-validation; not a feature)",
        "feature_groups": FEATURE_GROUPS,
        "descriptions": FEATURE_DESCRIPTIONS,
    })
    write_json(output_dir / "dataset_report.json", {
        "phase": 5,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "inputs": {p.relative_to(config.PROJECT_ROOT).as_posix(): sha256_file(p) for p in (
            cohort_path, config.PHASE2_DIR / "vep_annotations.csv",
            config.PHASE3_DIR / "variant_structure_features.csv", config.PHASE4_DIR / "foldx_ddg.csv")},
        "rows": len(rows),
        "genes": len({r["gene"] for r in rows}),
        "labels": dict(Counter(r["label"] for r in rows)),
        "features": len(ALL_FEATURES),
        "merge_losses": merge_losses,
        "missingness": missingness,
        "dataset_sha256": sha256_file(dataset_path),
    })
    write_json(output_dir / "validation_report.json", {"passed": passed, "checks": checks})

    print(f"{len(rows):,} rows, {len(ALL_FEATURES)} features")
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['check']}  {c['detail'] if not c['passed'] else ''}")
    return {"passed": passed}


def main(argv: list[str] | None = None) -> int:
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
