"""Download and validate only the five required AlphaFold structures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

COHORT = Path("data/processed/final_project_cohort.csv")
CANDIDATES = Path("data/processed/alphafold_structure_candidates.csv")
STRUCTURE_DIR = Path("data/raw/alphafold")
MAPPING = Path("data/processed/variant_structure_mapping.csv")
REPORT = Path("results/metrics/alphafold_structure_report.json")
GENES = {"BRCA1", "COL2A1", "GRIN2B", "NSD1", "TP53"}
AA_CODES = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
PROTEIN_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")


def _extract_cif_sequence(path: Path) -> str:
    data = MMCIF2Dict(str(path))
    for key in ("_entity_poly.pdbx_seq_one_letter_code_can", "_entity_poly.pdbx_seq_one_letter_code"):
        value = data.get(key)
        if value:
            sequence = value[0] if isinstance(value, list) else value
            return re.sub(r"\s+|\\n|\\r|\(|\)", "", sequence)
    raise ValueError(f"No polymer sequence found in CIF: {path}")


def _parse_variant(row: dict[str, str]) -> tuple[int, str, str] | None:
    match = PROTEIN_RE.search(row.get("Name", ""))
    if not match:
        return None
    return int(match.group(2)), AA_CODES[match.group(1)], AA_CODES[match.group(3)]


def acquire(
    cohort_path: Path = COHORT,
    candidates_path: Path = CANDIDATES,
    structure_dir: Path = STRUCTURE_DIR,
    mapping_path: Path = MAPPING,
    report_path: Path = REPORT,
) -> dict[str, Any]:
    candidates = {row["gene"]: row for row in csv.DictReader(candidates_path.open(encoding="utf-8", newline=""))}
    if not GENES.issubset(candidates) or any(candidates[gene]["structure_available"] != "True" for gene in GENES):
        raise ValueError("The AlphaFold candidates table does not contain five fully compatible genes")
    cohort_rows = list(csv.DictReader(cohort_path.open(encoding="utf-8", newline="")))
    for index, row in enumerate(cohort_rows, start=1):
        identity = "|".join([
            row["Assembly"], row["Chromosome"], row["Start"], row["Stop"],
            row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"], str(index),
        ])
        row["stable_variant_id"] = f"cohort_{index:05d}_{hashlib.sha1(identity.encode('utf-8')).hexdigest()[:12]}"
    structure_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "User-Agent": "BioDB-class-project/1.0"})
    gene_results: dict[str, dict[str, Any]] = {}
    mapping_rows: list[dict[str, Any]] = []

    for gene in sorted(GENES):
        accession = candidates[gene]["protein_accession"]
        metadata_response = session.get(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}", timeout=60)
        metadata_response.raise_for_status()
        metadata = metadata_response.json()[0]
        sequence = metadata["sequence"]
        structure_id = metadata["modelEntityId"]
        cif_url = metadata["cifUrl"]
        cif_path = structure_dir / f"{gene}_{accession}.cif"
        if not cif_path.exists() or cif_path.stat().st_size == 0:
            response = session.get(cif_url, timeout=180)
            response.raise_for_status()
            cif_path.write_bytes(response.content)
        cif_sequence = _extract_cif_sequence(cif_path)
        sequence_match = cif_sequence == sequence
        gene_variants = [row for row in cohort_rows if row["GeneSymbol"] == gene and row["Type"] != "Inversion"]
        compatible = 0
        out_of_range = 0
        sequence_mismatches = 0
        incompatible_details: list[dict[str, Any]] = []
        for row in gene_variants:
            parsed = _parse_variant(row)
            position = parsed[0] if parsed else None
            reference = parsed[1] if parsed else ""
            alternate = parsed[2] if parsed else ""
            reason = ""
            if not sequence_match:
                reason = "structure_sequence_mismatch_to_api_sequence"
                sequence_mismatches += 1
            elif parsed is None:
                reason = "missing_protein_notation"
            elif position > len(cif_sequence):
                reason = "out_of_range_position"
                out_of_range += 1
            elif cif_sequence[position - 1] != reference:
                reason = "reference_amino_acid_mismatch"
                sequence_mismatches += 1
            else:
                compatible += 1
            mapping_rows.append({
                "stable_variant_id": row.get("stable_variant_id", ""),
                "GeneSymbol": gene,
                "protein_change": re.search(r"\(p\.[^)]+\)", row["Name"]).group(0) if re.search(r"\(p\.[^)]+\)", row["Name"]) else "",
                "UniProt_accession": accession,
                "structure_id": structure_id,
                "residue_position": position or "",
                "reference_amino_acid": reference,
                "alternate_amino_acid": alternate,
                "structure_compatible": reason == "",
                "incompatibility_reason": reason,
            })
            if reason:
                incompatible_details.append({"stable_variant_id": row.get("stable_variant_id", ""), "reason": reason, "position": position})

        gene_results[gene] = {
            "gene": gene, "uniprot_accession": accession, "protein_length": len(sequence),
            "structure_id": structure_id, "structure_file": str(cif_path), "structure_format": "mmCIF",
            "download_url": cif_url, "source": "AlphaFold Protein Structure Database",
            "sequence_length_from_cif": len(cif_sequence), "sequence_matches_api": sequence_match,
            "variant_count": len(gene_variants), "compatible_variant_count": compatible,
            "coverage_percentage": round(compatible * 100 / len(gene_variants), 4) if gene_variants else 0.0,
            "out_of_range_count": out_of_range, "sequence_mismatch_count": sequence_mismatches,
            "incompatible_variants": incompatible_details,
            "isoform_rationale": "Reviewed human UniProt accession selected in Phase 4A; AlphaFold metadata identifies the same UniProt accession and sequence.",
            "file_sha256": hashlib.sha256(cif_path.read_bytes()).hexdigest(),
        }

    mapping_fields = ["stable_variant_id", "GeneSymbol", "protein_change", "UniProt_accession", "structure_id", "residue_position", "reference_amino_acid", "alternate_amino_acid", "structure_compatible", "incompatibility_reason"]
    with mapping_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=mapping_fields)
        writer.writeheader()
        writer.writerows(mapping_rows)

    report = {
        "access_date_utc": datetime.now(timezone.utc).isoformat(),
        "source": "AlphaFold Protein Structure Database",
        "source_api": "https://alphafold.ebi.ac.uk/api-docs",
        "cohort_file": str(cohort_path), "total_cohort_variants": len(cohort_rows),
        "downstream_variants_checked": len(mapping_rows), "inversion_variants_excluded": sum(row["Type"] == "Inversion" for row in cohort_rows),
        "genes": gene_results, "mapping_file": str(mapping_path),
        "structures_downloaded": sorted(str(path) for path in structure_dir.glob("*.cif")),
        "stability_calculations_performed": False,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Acquire and validate the five required AlphaFold structures.")
    parser.add_argument("--cohort", type=Path, default=COHORT)
    parser.add_argument("--candidates", type=Path, default=CANDIDATES)
    parser.add_argument("--structure-dir", type=Path, default=STRUCTURE_DIR)
    parser.add_argument("--mapping", type=Path, default=MAPPING)
    parser.add_argument("--report", type=Path, default=REPORT)
    args = parser.parse_args()
    report = acquire(args.cohort, args.candidates, args.structure_dir, args.mapping, args.report)
    print(f"Checked {report['downstream_variants_checked']:,} downstream variants across {len(report['genes'])} genes")
    for gene, result in report["genes"].items():
        print(gene, result["uniprot_accession"], result["structure_id"], result["coverage_percentage"])


if __name__ == "__main__":
    main()