"""Create the Phase 2 binary ClinVar modeling dataset.

The policy intentionally matches only the four exact aggregate classifications
observed in the downloaded ClinVar file. Composite, uncertain, conflicting,
missing, and non-clinical classifications are excluded.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

DEFAULT_INPUT = Path("data/raw/variant_summary.txt.gz")
DEFAULT_OUTPUT = Path("data/processed/clinvar_filtered.csv")
DEFAULT_REPORT = Path("results/metrics/clinvar_filtering_report.json")
MISSING_TOKENS = {"", "-", "NA", "N/A", "na", "n/a", "null", "None"}
LABELS = {
    "Pathogenic": 1,
    "Likely pathogenic": 1,
    "Benign": 0,
    "Likely benign": 0,
}
GENOMIC_COLUMNS = (
    "Assembly",
    "ChromosomeAccession",
    "Chromosome",
    "Start",
    "Stop",
)
PRIMARY_ALLELE_COLUMNS = ("ReferenceAllele", "AlternateAllele")
VCF_ALLELE_COLUMNS = ("ReferenceAlleleVCF", "AlternateAlleleVCF")


def _clean(value: str) -> str:
    return value.strip()


def _usable(value: str) -> bool:
    return _clean(value) not in MISSING_TOKENS


def _variant_alleles(row: dict[str, str]) -> tuple[str, str] | None:
    """Return VCF alleles when present, otherwise ClinVar primary alleles."""
    for reference_column, alternate_column in (
        VCF_ALLELE_COLUMNS,
        PRIMARY_ALLELE_COLUMNS,
    ):
        reference = _clean(row[reference_column])
        alternate = _clean(row[alternate_column])
        if _usable(reference) and _usable(alternate):
            return reference, alternate
    return None


def _variant_key(row: dict[str, str]) -> tuple[str, ...] | None:
    if not all(_usable(row[column]) for column in GENOMIC_COLUMNS):
        return None
    alleles = _variant_alleles(row)
    if alleles is None:
        return None
    return (
        _clean(row["Assembly"]),
        _clean(row["ChromosomeAccession"]),
        _clean(row["Chromosome"]),
        _clean(row["Start"]),
        _clean(row["Stop"]),
        *alleles,
    )


def _empty_report(input_path: Path, output_path: Path) -> dict[str, Any]:
    return {
        "source_file": str(input_path),
        "output_file": str(output_path),
        "classification_policy": {
            "pathogenicity_label_1": ["Pathogenic", "Likely pathogenic"],
            "pathogenicity_label_0": ["Benign", "Likely benign"],
            "excluded": "All other exact ClinicalSignificance values, including composite, uncertain, conflicting, missing, and non-clinical classifications.",
        },
        "variant_identifier": [
            "Assembly",
            "ChromosomeAccession",
            "Chromosome",
            "Start",
            "Stop",
            "ReferenceAlleleVCF or ReferenceAllele",
            "AlternateAlleleVCF or AlternateAllele",
        ],
        "rows_before_filtering": 0,
        "rows_after_filtering": 0,
        "pathogenic_count": 0,
        "benign_count": 0,
        "excluded_count": 0,
        "duplicate_count": 0,
        "class_balance": {"pathogenic": 0.0, "benign": 0.0},
        "excluded_by_reason": {
            "non_target_or_missing_classification": 0,
            "unusable_genomic_variant_information": 0,
            "duplicate_variant": 0,
        },
        "observed_target_classifications": sorted(LABELS),
    }


def process_clinvar(
    input_path: Path = DEFAULT_INPUT,
    output_path: Path = DEFAULT_OUTPUT,
    report_path: Path = DEFAULT_REPORT,
) -> dict[str, Any]:
    """Stream ClinVar, filter exact labels, deduplicate, and write a CSV."""
    if not input_path.exists():
        raise FileNotFoundError(f"ClinVar file not found: {input_path}")

    report = _empty_report(input_path, output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_suffix(output_path.suffix + ".tmp")
    database_path = output_path.with_suffix(".dedupe.sqlite3")
    temporary_output.unlink(missing_ok=True)
    database_path.unlink(missing_ok=True)
    classification_counts: Counter[str] = Counter()

    try:
        with gzip.open(input_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError("ClinVar file has no header")
            fieldnames = reader.fieldnames
            output_fields = [*fieldnames, "clinical_significance_original", "pathogenicity_label"]

            database = sqlite3.connect(database_path)
            try:
                database.execute("CREATE TABLE seen_variants (variant_key TEXT PRIMARY KEY)")
                with temporary_output.open("w", encoding="utf-8", newline="") as destination:
                    writer = csv.DictWriter(destination, fieldnames=output_fields)
                    writer.writeheader()

                    for row in reader:
                        report["rows_before_filtering"] += 1
                        classification = _clean(row.get("ClinicalSignificance", ""))
                        classification_counts[classification or "<MISSING>"] += 1
                        if classification not in LABELS:
                            report["excluded_by_reason"]["non_target_or_missing_classification"] += 1
                            continue

                        key = _variant_key(row)
                        if key is None:
                            report["excluded_by_reason"]["unusable_genomic_variant_information"] += 1
                            continue

                        key_text = "\x1f".join(key)
                        inserted = database.execute(
                            "INSERT OR IGNORE INTO seen_variants (variant_key) VALUES (?)",
                            (key_text,),
                        ).rowcount
                        if inserted == 0:
                            report["duplicate_count"] += 1
                            report["excluded_by_reason"]["duplicate_variant"] += 1
                            continue

                        output_row = dict(row)
                        output_row["clinical_significance_original"] = classification
                        output_row["pathogenicity_label"] = LABELS[classification]
                        writer.writerow(output_row)
                        report["rows_after_filtering"] += 1
                        if LABELS[classification] == 1:
                            report["pathogenic_count"] += 1
                        else:
                            report["benign_count"] += 1

                        if report["rows_after_filtering"] % 10000 == 0:
                            database.commit()
                database.commit()
            finally:
                database.close()

        temporary_output.replace(output_path)
    finally:
        temporary_output.unlink(missing_ok=True)
        database_path.unlink(missing_ok=True)

    report["excluded_count"] = report["rows_before_filtering"] - report["rows_after_filtering"]
    total = report["rows_after_filtering"]
    report["class_balance"] = {
        "pathogenic": round(report["pathogenic_count"] / total, 6) if total else 0.0,
        "benign": round(report["benign_count"] / total, 6) if total else 0.0,
    }
    report["classification_counts"] = dict(classification_counts)
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Create the binary ClinVar target dataset.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    report = process_clinvar(args.input, args.output, args.report)
    print(f"Rows before filtering: {report['rows_before_filtering']:,}")
    print(f"Rows after filtering: {report['rows_after_filtering']:,}")
    print(f"Pathogenic: {report['pathogenic_count']:,}")
    print(f"Benign: {report['benign_count']:,}")
    print(f"Duplicates removed: {report['duplicate_count']:,}")
    print(f"Filtered dataset written to: {args.output}")
    print(f"Report written to: {args.report}")


if __name__ == "__main__":
    main()