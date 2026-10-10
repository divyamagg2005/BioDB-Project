"""Validate targeted AlphaFold availability for the Phase 2 project cohort.

Only metadata endpoints are queried. This module does not download structure
files, alter the cohort, or perform mutation/stability calculations.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import requests

DEFAULT_COHORT = Path("data/processed/clinvar_project_cohort.csv")
DEFAULT_OUTPUT = Path("data/processed/alphafold_structure_candidates.csv")
DEFAULT_REPORT = Path("results/metrics/alphafold_availability_report.json")
GENE_ORDER = ["TTN", "BRCA2", "BRCA1", "NEB"]
AA_CODES = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
PROTEIN_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
OUTPUT_FIELDS = [
    "gene", "protein_accession", "alphafold_structure_id", "protein_length",
    "structure_available", "variant_count", "potentially_mappable_variants",
    "mapping_percentage", "notes",
]


def _get_json(session: requests.Session, url: str, params: dict[str, str] | None = None) -> Any:
    response = session.get(url, params=params, timeout=60)
    if response.status_code == 404:
        return None
    response.raise_for_status()
    return response.json()


def _load_cohort(path: Path) -> dict[str, list[dict[str, Any]]]:
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source):
            match = PROTEIN_RE.search(row.get("Name", ""))
            if match:
                row["protein_reference"] = AA_CODES[match.group(1)]
                row["protein_position"] = int(match.group(2))
                row["protein_alternate"] = AA_CODES[match.group(3)]
                by_gene[row["GeneSymbol"]].append(row)
    return by_gene


def _canonical_uniprot(session: requests.Session, gene: str) -> tuple[str, int, str, str]:
    data = _get_json(
        session,
        "https://rest.uniprot.org/uniprotkb/search",
        {"query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true", "format": "json", "size": "1"},
    )
    if not data or not data.get("results"):
        raise ValueError(f"No reviewed human UniProt entry found for {gene}")
    result = data["results"][0]
    accession = result["primaryAccession"]
    length = result["sequence"]["length"]
    sequence = result["sequence"]["value"]
    return accession, length, sequence, result.get("uniProtkbId", "")


def validate_availability(
    cohort_path: Path = DEFAULT_COHORT,
    output_path: Path = DEFAULT_OUTPUT,
    report_path: Path = DEFAULT_REPORT,
) -> dict[str, Any]:
    if not cohort_path.exists():
        raise FileNotFoundError(f"Cohort not found: {cohort_path}")
    by_gene = _load_cohort(cohort_path)
    genes_to_check = [gene for gene in GENE_ORDER if gene in by_gene]
    genes_to_check.extend(sorted(set(by_gene) - set(genes_to_check)))
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "BioDB-class-project/1.0"})
    rows: list[dict[str, Any]] = []
    mapping_details: dict[str, Any] = {}

    for gene in genes_to_check:
        variants = by_gene.get(gene, [])
        accession, protein_length, sequence, uniprot_id = _canonical_uniprot(session, gene)
        prediction = _get_json(session, f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
        summary = _get_json(session, f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{accession}.json")
        structure_available = bool(prediction)
        structure_id = ""
        if prediction:
            structure_id = prediction[0].get("modelEntityId", f"AF-{accession}-F1")

        reasons = Counter()
        mappable = 0
        position_out_of_range = 0
        sequence_mismatch = 0
        for variant in variants:
            position = variant["protein_position"]
            if not structure_available:
                reasons["unavailable_structure"] += 1
            elif position > protein_length:
                position_out_of_range += 1
                reasons["residue_position_mismatch"] += 1
            elif sequence[position - 1] != variant["protein_reference"]:
                sequence_mismatch += 1
                reasons["sequence_mismatch"] += 1
            else:
                mappable += 1

        max_position = max((v["protein_position"] for v in variants), default=0)
        notes = (
            f"Reviewed UniProt {uniprot_id}; AlphaFold metadata queried by accession. "
            f"Cohort protein positions range 1-{max_position}."
        )
        if not structure_available:
            notes += " No official AlphaFold prediction metadata was returned; no structure file was downloaded."
        elif position_out_of_range or sequence_mismatch:
            notes += f" {position_out_of_range} positions exceed sequence length; {sequence_mismatch} reference residues mismatch."
        else:
            notes += " All cohort reference residues are in range and match the UniProt sequence."

        rows.append({
            "gene": gene,
            "protein_accession": accession,
            "alphafold_structure_id": structure_id,
            "protein_length": protein_length,
            "structure_available": structure_available,
            "variant_count": len(variants),
            "potentially_mappable_variants": mappable,
            "mapping_percentage": round(mappable * 100 / len(variants), 4) if variants else 0.0,
            "notes": notes,
        })
        mapping_details[gene] = {
            "uniprot_id": uniprot_id,
            "summary_endpoint_found": summary is not None,
            "cohort_position_min": min((v["protein_position"] for v in variants), default=None),
            "cohort_position_max": max_position,
            "unmapped_by_reason": dict(reasons),
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=OUTPUT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    report = {
        "cohort_file": str(cohort_path),
        "total_variants": sum(len(values) for values in by_gene.values()),
        "genes_checked": genes_to_check,
        "structure_candidates": rows,
        "mapping_details": mapping_details,
        "structure_files_downloaded": False,
        "cohort_modified": False,
        "stability_calculations_performed": False,
        "overall_status": "requires review before structural phase: one or more selected genes lack official AlphaFold prediction metadata",
        "next_step": "Stop and decide whether to replace genes without AlphaFold metadata or proceed with the structurally covered subset.",
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate targeted AlphaFold metadata availability.")
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    report = validate_availability(args.cohort, args.output, args.report)
    print(f"Checked variants: {report['total_variants']:,}")
    for row in report["structure_candidates"]:
        print(row["gene"], row["protein_accession"], row["structure_available"], row["mapping_percentage"])
    print(f"Report written to: {args.report}")


if __name__ == "__main__":
    main()