"""Create the final AlphaFold-compatible ClinVar missense cohort."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

DEFAULT_INPUT = Path("data/processed/clinvar_filtered.csv")
DEFAULT_COMPATIBILITY = Path("data/processed/alphafold_structure_candidates.csv")
DEFAULT_OUTPUT = Path("data/processed/final_project_cohort.csv")
DEFAULT_GENES = Path("data/processed/final_selected_genes.csv")
DEFAULT_REPORT = Path("results/metrics/final_cohort_report.json")
FINAL_GENES = {"BRCA1", "COL2A1", "GRIN2B", "NSD1", "TP53"}
MISSING = {"", "-", "NA", "N/A", "na", "n/a", "null", "None"}
MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
AA_CODES = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
}
OUTPUT_FIELDS = [
    "#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
    "Assembly", "ChromosomeAccession", "Chromosome", "Start", "Stop",
    "ReferenceAlleleVCF", "AlternateAlleleVCF", "PositionVCF", "VariationID",
    "clinical_significance_original", "pathogenicity_label",
]


def _usable(value: str | None) -> bool:
    return (value or "").strip() not in MISSING


def _is_true_missense(row: dict[str, str]) -> bool:
    match = MISSENSE_RE.search(row.get("Name", ""))
    return (
        ":c." in row.get("Name", "")
        and match is not None
        and match.group(1) in AA_CODES
        and match.group(3) in AA_CODES
        and _usable(row.get("GeneSymbol"))
    )


def _sort_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row.get("GeneSymbol", ""), row.get("pathogenicity_label", ""),
        row.get("Assembly", ""), row.get("ChromosomeAccession", ""),
        row.get("Start", ""), row.get("Stop", ""),
        row.get("ReferenceAlleleVCF", ""), row.get("AlternateAlleleVCF", ""),
        row.get("#AlleleID", ""),
    )


def _load_compatibility(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as source:
        rows = {row["gene"]: row for row in csv.DictReader(source)}
    missing = FINAL_GENES - set(rows)
    if missing:
        raise ValueError(f"Missing AlphaFold compatibility records: {sorted(missing)}")
    invalid = [
        gene for gene in FINAL_GENES
        if rows[gene]["structure_available"] != "True"
        or float(rows[gene]["mapping_percentage"]) < 100.0
    ]
    if invalid:
        raise ValueError(f"Genes are not 100% AlphaFold-compatible: {sorted(invalid)}")
    return rows


def create_final_cohort(
    input_path: Path = DEFAULT_INPUT,
    compatibility_path: Path = DEFAULT_COMPATIBILITY,
    output_path: Path = DEFAULT_OUTPUT,
    genes_path: Path = DEFAULT_GENES,
    report_path: Path = DEFAULT_REPORT,
) -> dict[str, Any]:
    compatibility = _load_compatibility(compatibility_path)
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    counts = Counter()
    coding_count = 0
    missense_count = 0
    protein_count = 0
    before_count = 0

    with input_path.open("r", encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            before_count += 1
            if ":c." in row.get("Name", "") and _usable(row.get("GeneSymbol")):
                coding_count += 1
            if not _is_true_missense(row) or row["GeneSymbol"] not in FINAL_GENES:
                continue
            missense_count += 1
            protein_count += 1
            by_gene[row["GeneSymbol"]].append(row)
            counts[(row["GeneSymbol"], row["pathogenicity_label"])] += 1

    selected_rows: list[dict[str, str]] = []
    gene_rows: list[dict[str, Any]] = []
    excluded_genes = {
        "MECP2": "Poor AlphaFold mapping (9.9631%); excluded by final compatibility policy.",
        "COL7A1": "No official AlphaFold prediction metadata; excluded by final compatibility policy.",
        "DYNC1H1": "No official AlphaFold prediction metadata; excluded by final compatibility policy.",
        "PKD1": "No official AlphaFold prediction metadata; excluded by final compatibility policy.",
    }
    for gene in sorted(FINAL_GENES):
        rows = sorted(by_gene[gene], key=_sort_key)
        pathogenic_rows = [row for row in rows if row["pathogenicity_label"] == "1"]
        benign_rows = [row for row in rows if row["pathogenicity_label"] == "0"]
        per_class = min(len(pathogenic_rows), len(benign_rows))
        selected_rows.extend(pathogenic_rows[:per_class])
        selected_rows.extend(benign_rows[:per_class])
        gene_rows.append({
            "gene": gene,
            "protein_accession": compatibility[gene]["protein_accession"],
            "alphafold_structure_id": compatibility[gene]["alphafold_structure_id"],
            "available_pathogenic": len(pathogenic_rows),
            "available_benign": len(benign_rows),
            "selected_pathogenic": per_class,
            "selected_benign": per_class,
            "selected_total": 2 * per_class,
            "mapping_percentage": float(compatibility[gene]["mapping_percentage"]),
            "selection_policy": "Use the largest equal per-class count for each compatible gene after stable sorting; no random sampling.",
        })
    selected_rows.sort(key=_sort_key)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    genes_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_output.unlink(missing_ok=True)
    with temporary_output.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=OUTPUT_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected_rows)
    temporary_output.replace(output_path)

    with genes_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=list(gene_rows[0]))
        writer.writeheader()
        writer.writerows(gene_rows)

    final_counts = Counter(row["pathogenicity_label"] for row in selected_rows)
    report = {
        "source_file": str(input_path),
        "total_variants": len(selected_rows),
        "variants_before_filtering": before_count,
        "coding_variants_in_source": coding_count,
        "missense_candidates_in_selected_genes": missense_count,
        "pathogenic_count": final_counts["1"],
        "benign_count": final_counts["0"],
        "class_balance": {
            "pathogenic": round(final_counts["1"] / len(selected_rows), 6),
            "benign": round(final_counts["0"] / len(selected_rows), 6),
        },
        "variants_per_gene": {row["gene"]: row["selected_total"] for row in gene_rows},
        "variants_per_class_per_gene": {
            row["gene"]: {"pathogenic": row["selected_pathogenic"], "benign": row["selected_benign"]}
            for row in gene_rows
        },
        "protein_level_information_count": len(selected_rows),
        "alphafold_compatible_count": len(selected_rows),
        "selected_genes": sorted(FINAL_GENES),
        "excluded_genes_and_reasons": excluded_genes,
        "sampling_procedure": "For each of the five AlphaFold-compatible genes, sort by gene, label, assembly, genomic coordinates, alleles, and ClinVar allele ID; retain the same number from each class equal to that gene's smaller class count.",
        "target_note": "A 5,000-variant balanced cohort is unavailable from these five genes; this is the largest balanced cohort supported by their observed class counts.",
        "external_annotation_performed": False,
        "structures_downloaded": False,
        "stability_calculations_performed": False,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Create the final AlphaFold-compatible ClinVar cohort.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--compatibility", type=Path, default=DEFAULT_COMPATIBILITY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--genes", type=Path, default=DEFAULT_GENES)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    report = create_final_cohort(args.input, args.compatibility, args.output, args.genes, args.report)
    print(f"Final variants: {report['total_variants']:,}")
    print(f"Pathogenic: {report['pathogenic_count']:,}; benign: {report['benign_count']:,}")
    print(f"Final cohort written to: {args.output}")


if __name__ == "__main__":
    main()