"""Select a manageable, reproducible ClinVar missense research cohort.

The selector uses only the Phase 2 ClinVar CSV. It treats a row as a coding
missense candidate when its ClinVar name contains coding HGVS notation and a
three-letter amino-acid substitution such as ``(p.Arg330Met)``. Variants are
sorted by stable ClinVar identifiers before per-class selection; no random
sampling or external annotation is used.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

DEFAULT_INPUT = Path("data/processed/clinvar_filtered.csv")
DEFAULT_OUTPUT = Path("data/processed/clinvar_project_cohort.csv")
DEFAULT_GENES = Path("data/processed/selected_genes.csv")
DEFAULT_REPORT = Path("results/metrics/clinvar_project_cohort_report.json")
GENE_LIMIT = 10
PER_CLASS_PER_GENE = 750
MIN_CLASS_SUPPORT = 25
MISSING = {"", "-", "NA", "N/A", "na", "n/a", "null", "None"}
MISSENSE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
AA_CODES = {
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
}


def _usable(value: str | None) -> bool:
    return (value or "").strip() not in MISSING


def _is_missense(row: dict[str, str]) -> bool:
    name = row.get("Name", "")
    match = MISSENSE_RE.search(name)
    return (
        ":c." in name
        and match is not None
        and match.group(1) in AA_CODES
        and match.group(3) in AA_CODES
    )


def _is_coding(row: dict[str, str]) -> bool:
    return ":c." in row.get("Name", "") and _usable(row.get("GeneSymbol"))


def _variant_sort_key(row: dict[str, str]) -> tuple[str, ...]:
    return (
        row.get("GeneSymbol", ""),
        row.get("pathogenicity_label", ""),
        row.get("Assembly", ""),
        row.get("ChromosomeAccession", ""),
        row.get("Start", ""),
        row.get("Stop", ""),
        row.get("ReferenceAlleleVCF", ""),
        row.get("AlternateAlleleVCF", ""),
        row.get("#AlleleID", ""),
    )


def _scan_candidates(input_path: Path) -> tuple[dict[str, Any], dict[str, list[dict[str, str]]]]:
    summary: dict[str, Any] = {
        "variants_before_filtering": 0,
        "variants_after_clinical_label_filtering": 0,
        "coding_variants": 0,
        "missense_variants": 0,
        "protein_level_variants": 0,
        "excluded_without_protein_mapping": 0,
    }
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    with input_path.open("r", encoding="utf-8", newline="") as source:
        reader = csv.DictReader(source)
        for row in reader:
            summary["variants_before_filtering"] += 1
            summary["variants_after_clinical_label_filtering"] += 1
            if not _is_coding(row):
                continue
            summary["coding_variants"] += 1
            if not _is_missense(row):
                continue
            summary["missense_variants"] += 1
            summary["protein_level_variants"] += 1
            by_gene[row["GeneSymbol"].strip()].append(row)
    summary["candidate_gene_count"] = len(by_gene)
    return summary, by_gene


def _rank_genes(by_gene: dict[str, list[dict[str, str]]]) -> list[dict[str, Any]]:
    ranked = []
    for gene, rows in by_gene.items():
        counts = Counter(row["pathogenicity_label"] for row in rows)
        pathogenic = counts.get("1", 0)
        benign = counts.get("0", 0)
        if pathogenic < MIN_CLASS_SUPPORT or benign < MIN_CLASS_SUPPORT:
            continue
        minimum = min(pathogenic, benign)
        balance = minimum / max(pathogenic, benign)
        ranked.append({
            "gene": gene,
            "pathogenic_count": pathogenic,
            "benign_count": benign,
            "total_count": pathogenic + benign,
            "minority_class_count": minimum,
            "class_balance": round(balance, 6),
            "selection_policy": "Rank by minority class support, then class balance, then total support, then gene symbol.",
        })
    return sorted(
        ranked,
        key=lambda item: (-item["minority_class_count"], -item["class_balance"], -item["total_count"], item["gene"]),
    )


def _select_genes(ranked: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    projected = 0
    for candidate in ranked[:GENE_LIMIT]:
        selected.append(candidate)
        projected += 2 * min(PER_CLASS_PER_GENE, candidate["minority_class_count"])
        if projected >= 5000:
            break
    return selected


def create_cohort(
    input_path: Path = DEFAULT_INPUT,
    output_path: Path = DEFAULT_OUTPUT,
    genes_path: Path = DEFAULT_GENES,
    report_path: Path = DEFAULT_REPORT,
) -> dict[str, Any]:
    if not input_path.exists():
        raise FileNotFoundError(f"Phase 2 dataset not found: {input_path}")
    summary, by_gene = _scan_candidates(input_path)
    ranked = _rank_genes(by_gene)
    selected = _select_genes(ranked)
    selected_gene_names = {item["gene"] for item in selected}
    if not selected:
        raise ValueError("No genes met the minimum pathogenic and benign support threshold")

    selected_rows: list[dict[str, str]] = []
    for gene in selected_gene_names:
        rows = sorted(by_gene[gene], key=_variant_sort_key)
        for label in ("0", "1"):
            class_rows = [row for row in rows if row["pathogenicity_label"] == label]
            selected_rows.extend(class_rows[:PER_CLASS_PER_GENE])
    selected_rows.sort(key=_variant_sort_key)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    genes_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_output = output_path.with_suffix(output_path.suffix + ".tmp")
    temporary_output.unlink(missing_ok=True)
    with output_path.parent.joinpath(".cohort_genes.tmp").open("w", encoding="utf-8") as marker:
        marker.write("selection in progress")

    output_fields = [
        "#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "ClinicalSignificance",
        "Assembly", "ChromosomeAccession", "Chromosome", "Start", "Stop",
        "ReferenceAlleleVCF", "AlternateAlleleVCF", "PositionVCF", "VariationID",
        "clinical_significance_original", "pathogenicity_label",
    ]
    with temporary_output.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=output_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(selected_rows)
    temporary_output.replace(output_path)
    output_path.parent.joinpath(".cohort_genes.tmp").unlink(missing_ok=True)

    final_counts = Counter(row["pathogenicity_label"] for row in selected_rows)
    with genes_path.open("w", encoding="utf-8", newline="") as destination:
        fields = ["gene", "pathogenic_count", "benign_count", "total_count", "minority_class_count", "class_balance", "selection_policy"]
        writer = csv.DictWriter(destination, fieldnames=fields)
        writer.writeheader()
        writer.writerows(selected)

    summary.update({
        "candidate_genes": ranked,
        "final_selected_genes": selected,
        "final_variant_count": len(selected_rows),
        "final_pathogenic_count": final_counts.get("1", 0),
        "final_benign_count": final_counts.get("0", 0),
        "final_class_balance": {
            "pathogenic": round(final_counts.get("1", 0) / len(selected_rows), 6),
            "benign": round(final_counts.get("0", 0) / len(selected_rows), 6),
        },
        "protein_level_information_count": len(selected_rows),
        "sampling_procedure": {
            "gene_limit": GENE_LIMIT,
            "per_class_per_gene_limit": PER_CLASS_PER_GENE,
            "minimum_class_support_for_candidate": MIN_CLASS_SUPPORT,
            "variant_order": "Stable lexical order of gene, label, assembly, genomic coordinates, alleles, and ClinVar allele ID.",
            "random_sampling": False,
        },
        "alphafold_status": "Not queried in Phase 2 revision; selected genes have usable protein-level ClinVar notation and are candidates for targeted AlphaFold lookup in a later phase.",
    })
    report_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Create the manageable ClinVar protein-variant cohort.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--genes", type=Path, default=DEFAULT_GENES)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    report = create_cohort(args.input, args.output, args.genes, args.report)
    print(f"Missense variants: {report['missense_variants']:,}")
    print(f"Selected genes: {', '.join(item['gene'] for item in report['final_selected_genes'])}")
    print(f"Final cohort: {report['final_variant_count']:,}")
    print(f"Cohort written to: {args.output}")
    print(f"Gene selection written to: {args.genes}")
    print(f"Report written to: {args.report}")


if __name__ == "__main__":
    main()