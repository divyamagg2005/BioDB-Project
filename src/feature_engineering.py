"""Prepare a fixed, ClinVar-only variant feature table for later ML work.

This phase does not download annotation databases and does not invent values
for unavailable prediction scores. Unavailable score columns are retained as
missing values so the feature contract is explicit and reproducible.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

DEFAULT_INPUT = Path("data/processed/clinvar_filtered.csv")
DEFAULT_OUTPUT = Path("data/processed/variant_features.csv")
DEFAULT_REPORT = Path("results/metrics/variant_features_report.json")
MISSING_TOKENS = {"", "-", "NA", "N/A", "na", "n/a", "null", "None"}

FEATURE_COLUMNS = [
    "clinvar_allele_id",
    "assembly",
    "chromosome",
    "chromosome_accession",
    "start",
    "stop",
    "position_vcf",
    "reference_allele",
    "alternate_allele",
    "variant_type",
    "gene_id",
    "gene_symbol",
    "cadd_score",
    "sift_score",
    "polyphen2_score",
    "revel_score",
    "alphamissense_score",
    "allele_frequency",
    "conservation_score",
    "clinical_significance_original",
    "pathogenicity_label",
]

FEATURE_METADATA: dict[str, dict[str, str]] = {
    "clinvar_allele_id": {"source": "ClinVar #AlleleID", "meaning": "ClinVar allele identifier", "data_type": "integer-like identifier", "preprocessing": "Treat as an identifier, not a continuous measurement."},
    "assembly": {"source": "ClinVar Assembly", "meaning": "Reference genome assembly", "data_type": "categorical", "preprocessing": "Encode categorically or split by assembly."},
    "chromosome": {"source": "ClinVar Chromosome", "meaning": "Chromosome label", "data_type": "categorical", "preprocessing": "Normalize labels and encode categorically."},
    "chromosome_accession": {"source": "ClinVar ChromosomeAccession", "meaning": "Reference sequence accession", "data_type": "categorical", "preprocessing": "Encode categorically if used."},
    "start": {"source": "ClinVar Start", "meaning": "1-based genomic start coordinate", "data_type": "integer", "preprocessing": "Convert to numeric; verify assembly before modeling."},
    "stop": {"source": "ClinVar Stop", "meaning": "1-based genomic stop coordinate", "data_type": "integer", "preprocessing": "Convert to numeric; verify assembly before modeling."},
    "position_vcf": {"source": "ClinVar PositionVCF", "meaning": "VCF-style genomic position", "data_type": "integer", "preprocessing": "Convert to numeric and retain coordinate convention explicitly."},
    "reference_allele": {"source": "ClinVar ReferenceAlleleVCF, fallback ReferenceAllele", "meaning": "Reference sequence allele", "data_type": "string", "preprocessing": "Use VCF allele when available; encode sequence length/content later."},
    "alternate_allele": {"source": "ClinVar AlternateAlleleVCF, fallback AlternateAllele", "meaning": "Alternate sequence allele", "data_type": "string", "preprocessing": "Use VCF allele when available; encode sequence length/content later."},
    "variant_type": {"source": "ClinVar Type", "meaning": "ClinVar variant category", "data_type": "categorical", "preprocessing": "Normalize and one-hot encode."},
    "gene_id": {"source": "ClinVar GeneID", "meaning": "NCBI gene identifier", "data_type": "integer-like identifier", "preprocessing": "Treat as categorical or identifier."},
    "gene_symbol": {"source": "ClinVar GeneSymbol", "meaning": "Associated gene symbol", "data_type": "categorical", "preprocessing": "Normalize and encode categorically."},
    "cadd_score": {"source": "Unavailable in current ClinVar data", "meaning": "CADD deleteriousness score", "data_type": "float", "preprocessing": "Remain missing; obtain only from a documented lightweight annotation source later."},
    "sift_score": {"source": "Unavailable in current ClinVar data", "meaning": "SIFT predicted functional impact", "data_type": "float", "preprocessing": "Remain missing; do not impute as a biological value."},
    "polyphen2_score": {"source": "Unavailable in current ClinVar data", "meaning": "PolyPhen-2 predicted functional impact", "data_type": "float", "preprocessing": "Remain missing; do not impute as a biological value."},
    "revel_score": {"source": "Unavailable in current ClinVar data", "meaning": "REVEL ensemble missense score", "data_type": "float", "preprocessing": "Remain missing; do not impute as a biological value."},
    "alphamissense_score": {"source": "Unavailable in current ClinVar data", "meaning": "AlphaMissense pathogenicity score", "data_type": "float", "preprocessing": "Remain missing; no AlphaFold/AlphaMissense download in this phase."},
    "allele_frequency": {"source": "Unavailable in current ClinVar data", "meaning": "Population allele frequency", "data_type": "float", "preprocessing": "Remain missing; no population database was added."},
    "conservation_score": {"source": "Unavailable in current ClinVar data", "meaning": "Evolutionary conservation measure", "data_type": "float", "preprocessing": "Remain missing; no conservation database was added."},
    "clinical_significance_original": {"source": "Phase 2 ClinVar clinical significance", "meaning": "Original retained clinical classification", "data_type": "categorical", "preprocessing": "Keep for provenance; never use as an input feature without leakage review."},
    "pathogenicity_label": {"source": "Phase 2 exact-label policy", "meaning": "Binary target: pathogenicity", "data_type": "integer 0/1", "preprocessing": "Target column, not an input feature."},
}


def _clean(value: str | None) -> str:
    return (value or "").strip()


def _usable(value: str | None) -> bool:
    return _clean(value) not in MISSING_TOKENS


def _first_usable(row: dict[str, str], *columns: str) -> str:
    for column in columns:
        if _usable(row.get(column)):
            return _clean(row[column])
    return ""


def _feature_row(row: dict[str, str]) -> dict[str, str]:
    return {
        "clinvar_allele_id": _clean(row.get("#AlleleID")),
        "assembly": _clean(row.get("Assembly")),
        "chromosome": _clean(row.get("Chromosome")),
        "chromosome_accession": _clean(row.get("ChromosomeAccession")),
        "start": _clean(row.get("Start")),
        "stop": _clean(row.get("Stop")),
        "position_vcf": _clean(row.get("PositionVCF")),
        "reference_allele": _first_usable(row, "ReferenceAlleleVCF", "ReferenceAllele"),
        "alternate_allele": _first_usable(row, "AlternateAlleleVCF", "AlternateAllele"),
        "variant_type": _clean(row.get("Type")),
        "gene_id": _clean(row.get("GeneID")),
        "gene_symbol": _clean(row.get("GeneSymbol")),
        "cadd_score": "",
        "sift_score": "",
        "polyphen2_score": "",
        "revel_score": "",
        "alphamissense_score": "",
        "allele_frequency": "",
        "conservation_score": "",
        "clinical_significance_original": _clean(row.get("clinical_significance_original")),
        "pathogenicity_label": _clean(row.get("pathogenicity_label")),
    }


def generate_features(
    input_path: Path = DEFAULT_INPUT,
    output_path: Path = DEFAULT_OUTPUT,
    report_path: Path = DEFAULT_REPORT,
) -> dict[str, Any]:
    """Stream the Phase 2 dataset into the fixed Phase 3 feature schema."""
    if not input_path.exists():
        raise FileNotFoundError(f"Filtered ClinVar file not found: {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_output.unlink(missing_ok=True)
    row_count = 0
    missing_counts = Counter()
    label_counts = Counter()

    with input_path.open("r", encoding="utf-8", newline="") as source, temporary_output.open("w", encoding="utf-8", newline="") as destination:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError("Feature input has no header")
        missing_input_columns = {"#AlleleID", "Assembly", "Chromosome", "Start", "Stop", "Type", "GeneSymbol", "pathogenicity_label"} - set(reader.fieldnames)
        if missing_input_columns:
            raise ValueError(f"Missing required Phase 2 columns: {sorted(missing_input_columns)}")
        writer = csv.DictWriter(destination, fieldnames=FEATURE_COLUMNS)
        writer.writeheader()
        for row in reader:
            features = _feature_row(row)
            writer.writerow(features)
            row_count += 1
            label_counts[features["pathogenicity_label"]] += 1
            for column, value in features.items():
                if not _usable(value):
                    missing_counts[column] += 1

    temporary_output.replace(output_path)
    missingness = {
        column: {
            "count": missing_counts[column],
            "percentage": round(missing_counts[column] * 100 / row_count, 4) if row_count else 0.0,
        }
        for column in FEATURE_COLUMNS
    }
    report = {
        "source_file": str(input_path),
        "output_file": str(output_path),
        "row_count": row_count,
        "feature_count": len(FEATURE_COLUMNS),
        "features": FEATURE_METADATA,
        "missingness": missingness,
        "label_counts": dict(label_counts),
        "unavailable_requested_features": [
            "cadd_score", "sift_score", "polyphen2_score", "revel_score",
            "alphamissense_score", "allele_frequency", "conservation_score",
        ],
        "source_policy": "ClinVar-only for Phase 3; no dbNSFP, AlphaMissense, AlphaFold, or external annotation download was performed.",
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Create the fixed ClinVar variant feature table.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    report = generate_features(args.input, args.output, args.report)
    print(f"Feature rows: {report['row_count']:,}")
    print(f"Features: {report['feature_count']}")
    print(f"Feature table written to: {args.output}")
    print(f"Report written to: {args.report}")


if __name__ == "__main__":
    main()