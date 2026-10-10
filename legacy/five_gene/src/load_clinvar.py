"""Inspect the downloaded ClinVar variant summary in chunks.

This module does not alter the raw ClinVar archive or perform final label
filtering. It produces a compact JSON report with schema and distributions.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

DEFAULT_INPUT = Path("data/raw/variant_summary.txt.gz")
DEFAULT_OUTPUT = Path("results/metrics/clinvar_report.json")
DEFAULT_CHUNK_SIZE = 100_000
MISSING_TOKENS = {"", "-", "NA", "N/A", "na", "n/a", "null", "None"}


def _is_missing(value: str) -> bool:
    return value.strip() in MISSING_TOKENS


def _observed_type(value: str) -> str:
    """Classify a non-missing field value without coercing its contents."""
    value = value.strip()
    try:
        int(value)
        return "int64"
    except ValueError:
        try:
            float(value)
            return "float64"
        except ValueError:
            return "object"


def _add_distribution(target: dict[str, int], value: str) -> None:
    key = "<MISSING>" if _is_missing(value) else value
    target[key] = target.get(key, 0) + 1


def generate_report(
    input_path: Path = DEFAULT_INPUT,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, Any]:
    """Read ClinVar in chunks and return a compact dataset report."""
    if not input_path.exists():
        raise FileNotFoundError(f"ClinVar file not found: {input_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero")

    with gzip.open(input_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        columns = next(reader)

    missing_counts = {column: 0 for column in columns}
    observed_dtypes: dict[str, set[str]] = defaultdict(set)
    column_indexes = {column: index for index, column in enumerate(columns)}
    distributions = {
        "ClinicalSignificance": {},
        "Type": {},
        "Chromosome": {},
        "ReviewStatus": {},
    }
    row_count = 0
    malformed_row_count = 0

    with gzip.open(input_path, "rt", encoding="utf-8", errors="replace", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        next(reader)
        for row in reader:
            row_count += 1
            if len(row) != len(columns):
                malformed_row_count += 1
            values = row[: len(columns)]
            if len(values) < len(columns):
                values.extend([""] * (len(columns) - len(values)))
            for index, column in enumerate(columns):
                value = values[index]
                if _is_missing(value):
                    missing_counts[column] += 1
                else:
                    observed_dtypes[column].add(_observed_type(value))
            for column, distribution in distributions.items():
                if column in column_indexes:
                    _add_distribution(distribution, values[column_indexes[column]])

    file_size_bytes = input_path.stat().st_size
    data_types = {
        column: {
            "inferred": (
                next(iter(dtypes)) if len(dtypes) == 1 else "mixed"
            ),
            "observed": sorted(dtypes),
        }
        for column, dtypes in observed_dtypes.items()
    }
    missing_values = {
        column: {
            "count": count,
            "percentage": round((count / row_count) * 100, 4)
            if row_count
            else 0.0,
        }
        for column, count in missing_counts.items()
    }

    return {
        "source_file": str(input_path),
        "file_size_bytes": file_size_bytes,
        "chunk_size": chunk_size,
        "row_count": row_count,
        "column_count": len(columns),
        "columns": columns,
        "malformed_row_count": malformed_row_count,
        "data_types": data_types,
        "clinical_significance_distribution": distributions["ClinicalSignificance"],
        "variant_type_distribution": distributions["Type"],
        "chromosome_distribution": distributions["Chromosome"],
        "review_status_distribution": distributions["ReviewStatus"],
        "missing_values": missing_values,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a chunked inspection report for a ClinVar archive."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Path to the raw ClinVar .txt.gz file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Path for the generated JSON report.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help="Number of rows processed per chunk.",
    )
    args = parser.parse_args()

    report = generate_report(args.input, args.chunk_size)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"ClinVar rows: {report['row_count']:,}")
    print(f"ClinVar columns: {report['column_count']}")
    print(f"Report written to: {args.output}")


if __name__ == "__main__":
    main()