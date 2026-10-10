"""Phase 3: AlphaFold Structure Integration for 100-gene cohort.

This module:
- Downloads AlphaFold structures for 100 genes
- Maps variants to protein residues
- Validates UniProt → AlphaFold → residue consistency
- Creates variant-structure mapping table

Requirements:
- Input: Phase 1 gene list (phase1_selected_genes.csv)
- Input: Phase 2 annotated cohort (phase2_annotated_cohort.csv)
- Output: AlphaFold structures + variant mapping table

Source: AlphaFold Protein Structure Database (NOT PDB)
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_GENES = Path("data/processed/phase1/phase1_selected_genes.csv")
DEFAULT_COHORT = Path("data/processed/phase2/phase2_annotated_cohort.csv")
DEFAULT_STRUCTURE_DIR = Path("data/raw/alphafold_phase3")
DEFAULT_MAPPING = Path("data/processed/phase3/phase3_variant_structure_mapping.csv")
DEFAULT_REPORT = Path("data/processed/phase3/phase3_alphafold_report.json")
DEFAULT_VALIDATION = Path("data/processed/phase3/phase3_validation_report.json")

# API endpoints
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"

# Rate limiting
REQUEST_DELAY = 0.5  # seconds between API calls
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds

# Amino acid codes
AA_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# Protein change regex (p.XXXnnnYYY)
PROTEIN_CHANGE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def parse_protein_change(name: str) -> tuple[int, str, str] | None:
    """Extract position and amino acids from protein change notation.
    
    Args:
        name: Variant name like "NM_000033.4(ABCD1):c.19C>T (p.Pro7Ser)"
    
    Returns:
        (position, reference_aa, alternate_aa) or None if not parseable
    """
    match = PROTEIN_CHANGE_RE.search(name)
    if not match:
        return None
    
    ref_aa_name = match.group(1)
    position = int(match.group(2))
    alt_aa_name = match.group(3)
    
    ref_aa = AA_THREE_TO_ONE.get(ref_aa_name)
    alt_aa = AA_THREE_TO_ONE.get(alt_aa_name)
    
    if ref_aa is None or alt_aa is None:
        return None
    
    return position, ref_aa, alt_aa


def extract_cif_sequence(cif_path: Path) -> str:
    """Extract one-letter sequence from mmCIF file.
    
    Args:
        cif_path: Path to AlphaFold .cif file
    
    Returns:
        One-letter amino acid sequence
    """
    data = MMCIF2Dict(str(cif_path))
    
    # Try canonical sequence first
    for key in ("_entity_poly.pdbx_seq_one_letter_code_can", "_entity_poly.pdbx_seq_one_letter_code"):
        value = data.get(key)
        if value:
            sequence = value[0] if isinstance(value, list) else value
            # Remove whitespace and line continuation markers
            return re.sub(r"\s+|\\n|\\r|\(|\)", "", sequence)
    
    raise ValueError(f"No polymer sequence found in CIF: {cif_path}")


def download_alphafold_structure(
    session: requests.Session,
    uniprot_id: str,
    gene: str,
    structure_dir: Path,
    gene_index: int,
    total_genes: int,
) -> dict[str, Any]:
    """Download and validate AlphaFold structure for one gene.
    
    Args:
        session: Requests session
        uniprot_id: UniProt accession (e.g., "P38398")
        gene: Gene symbol
        structure_dir: Directory to save structures
        gene_index: Current gene index (for progress)
        total_genes: Total number of genes
    
    Returns:
        Dict with structure metadata and validation results
    """
    cif_path = structure_dir / f"{gene}_{uniprot_id}.cif"
    
    # Get metadata from AlphaFold API
    api_url = ALPHAFOLD_API.format(uniprot_id=uniprot_id)
    
    for attempt in range(1, MAX_RETRIES + 1):
        if attempt > 1:
            time.sleep(RETRY_DELAY)
        
        try:
            print(f"  [{gene_index}/{total_genes}] {gene} ({uniprot_id}) - attempt {attempt}/{MAX_RETRIES}")
            
            # Get prediction metadata
            metadata_response = session.get(api_url, timeout=60)
            
            if metadata_response.status_code == 404:
                return {
                    "gene": gene,
                    "uniprot_id": uniprot_id,
                    "status": "not_found",
                    "error": "AlphaFold structure not available for this UniProt ID",
                }
            
            metadata_response.raise_for_status()
            predictions = metadata_response.json()
            
            if not predictions:
                return {
                    "gene": gene,
                    "uniprot_id": uniprot_id,
                    "status": "no_predictions",
                    "error": "Empty response from AlphaFold API",
                }
            
            # Take first prediction (reviewed canonical isoform)
            prediction = predictions[0]
            sequence = prediction["sequence"]
            structure_id = prediction["modelEntityId"]
            cif_url = prediction["cifUrl"]
            
            # Download structure if not present
            if not cif_path.exists() or cif_path.stat().st_size == 0:
                print(f"    Downloading {cif_url}")
                cif_response = session.get(cif_url, timeout=180)
                cif_response.raise_for_status()
                cif_path.write_bytes(cif_response.content)
                time.sleep(REQUEST_DELAY)
            
            # Validate structure
            cif_sequence = extract_cif_sequence(cif_path)
            sequence_match = cif_sequence == sequence
            
            # Compute SHA256
            file_sha256 = hashlib.sha256(cif_path.read_bytes()).hexdigest()
            
            return {
                "gene": gene,
                "uniprot_id": uniprot_id,
                "status": "success",
                "structure_id": structure_id,
                "structure_file": str(cif_path),
                "structure_format": "mmCIF",
                "download_url": cif_url,
                "protein_length_api": len(sequence),
                "protein_length_cif": len(cif_sequence),
                "sequence_match": sequence_match,
                "file_sha256": file_sha256,
                "source": "AlphaFold Protein Structure Database",
            }
            
        except requests.RequestException as e:
            error_msg = str(e)[:200]
            if attempt == MAX_RETRIES:
                return {
                    "gene": gene,
                    "uniprot_id": uniprot_id,
                    "status": "download_failed",
                    "error": error_msg,
                }
        except Exception as e:
            return {
                "gene": gene,
                "uniprot_id": uniprot_id,
                "status": "error",
                "error": str(e)[:200],
            }
    
    return {
        "gene": gene,
        "uniprot_id": uniprot_id,
        "status": "max_retries_exceeded",
        "error": "Failed after maximum retries",
    }


def validate_variant_mapping(
    variant: dict[str, str],
    structure_data: dict[str, Any],
    cached_sequence: str = "",
) -> tuple[bool, str, int | None, str, str]:
    """Validate a single variant against structure data.
    
    Args:
        variant: Variant row from cohort
        structure_data: Structure metadata dict
        cached_sequence: Pre-parsed sequence from CIF (avoids re-parsing)
    
    Returns:
        (is_valid, reason, position, ref_aa, alt_aa)
    """
    # Check structure was downloaded
    if structure_data.get("status") != "success":
        return False, f"structure_unavailable:{structure_data.get('status')}", None, "", ""
    
    # Parse protein change
    parsed = parse_protein_change(variant.get("Name", ""))
    if parsed is None:
        return False, "missing_protein_notation", None, "", ""
    
    position, ref_aa, alt_aa = parsed
    
    # Check sequence match
    if not structure_data.get("sequence_match", False):
        return False, "structure_sequence_mismatch", position, ref_aa, alt_aa
    
    # Use cached sequence or load from CIF
    if cached_sequence:
        cif_sequence = cached_sequence
    else:
        cif_path = Path(structure_data["structure_file"])
        try:
            cif_sequence = extract_cif_sequence(cif_path)
        except Exception:
            return False, "cif_parse_error", position, ref_aa, alt_aa
    
    # Check position in range
    if position > len(cif_sequence):
        return False, f"out_of_range:pos={position}>len={len(cif_sequence)}", position, ref_aa, alt_aa
    
    # Check reference amino acid
    actual_ref = cif_sequence[position - 1]  # 1-indexed position
    if actual_ref != ref_aa:
        return False, f"reference_mismatch:expected={ref_aa},actual={actual_ref}", position, ref_aa, alt_aa
    
    return True, "valid", position, ref_aa, alt_aa


# ============================================================================
# MAIN INTEGRATION FUNCTION
# ============================================================================

def integrate_alphafold(
    genes_path: Path = DEFAULT_GENES,
    cohort_path: Path = DEFAULT_COHORT,
    structure_dir: Path = DEFAULT_STRUCTURE_DIR,
    mapping_path: Path = DEFAULT_MAPPING,
    report_path: Path = DEFAULT_REPORT,
    validation_path: Path = DEFAULT_VALIDATION,
) -> dict[str, Any]:
    """Run Phase 3 AlphaFold integration.
    
    Returns:
        Report dictionary with all statistics
    """
    print("=" * 80)
    print("PHASE 3: ALPHAFOLD STRUCTURE INTEGRATION")
    print("=" * 80)
    
    # Create directories
    structure_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Load genes
    print("\nLoading gene list...")
    genes_data: list[dict[str, str]] = []
    with genes_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        genes_data = list(reader)
    
    print(f"  Genes: {len(genes_data)}")
    print(f"  UniProt IDs: {', '.join(g['uniprot_accession'] for g in genes_data[:5])}...")
    
    # Load cohort
    print("\nLoading annotated cohort...")
    cohort_data: list[dict[str, str]] = []
    with cohort_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        cohort_data = list(reader)
    
    print(f"  Variants: {len(cohort_data)}")
    
    # Group variants by gene
    variants_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in cohort_data:
        gene = row.get("GeneSymbol", "")
        if gene:
            variants_by_gene[gene].append(row)
    
    print(f"  Genes in cohort: {len(variants_by_gene)}")
    
    # Setup session
    session = requests.Session()
    session.headers.update({
        "Accept": "application/json",
        "User-Agent": "BioDB-Phase3-AlphaFold/1.0"
    })
    
    # Download structures
    print("\nDownloading AlphaFold structures...")
    structure_results: dict[str, dict[str, Any]] = {}
    
    for idx, gene_row in enumerate(genes_data, start=1):
        gene = gene_row["gene"]
        uniprot_id = gene_row["uniprot_accession"]
        
        result = download_alphafold_structure(
            session, uniprot_id, gene, structure_dir, idx, len(genes_data)
        )
        structure_results[gene] = result
        
        # Progress update
        if result["status"] == "success":
            print(f"    ✓ {gene}: {result['protein_length_cif']} aa")
        else:
            print(f"    ✗ {gene}: {result['status']}")
        
        time.sleep(REQUEST_DELAY)
    
    # Count successes/failures
    successful = sum(1 for r in structure_results.values() if r["status"] == "success")
    failed = sum(1 for r in structure_results.values() if r["status"] != "success")
    
    print(f"\nStructure download results:")
    print(f"  Successful: {successful}")
    print(f"  Failed: {failed}")
    
    # Map variants (with sequence caching for performance)
    print("\nMapping variants to structures...")
    mapping_rows: list[dict[str, Any]] = []
    
    valid_count = 0
    invalid_reasons: Counter = Counter()
    
    # Pre-load sequences for all genes (cache to avoid re-parsing)
    print("  Caching sequences...")
    sequence_cache: dict[str, str] = {}
    for gene, structure_data in structure_results.items():
        if structure_data.get("status") == "success":
            cif_path = Path(structure_data["structure_file"])
            try:
                sequence_cache[gene] = extract_cif_sequence(cif_path)
            except Exception:
                sequence_cache[gene] = ""
    
    print(f"  Cached {len(sequence_cache)} sequences")
    print("  Mapping variants...")
    
    for gene, variants in variants_by_gene.items():
        structure_data = structure_results.get(gene, {"status": "no_gene_in_list"})
        cached_sequence = sequence_cache.get(gene, "")
        
        for idx, variant in enumerate(variants):
            if idx % 100 == 0:
                print(f"    {gene}: {idx}/{len(variants)} variants")
            
            is_valid, reason, position, ref_aa, alt_aa = validate_variant_mapping(
                variant, structure_data, cached_sequence
            )
            
            if is_valid:
                valid_count += 1
            else:
                invalid_reasons[reason.split(":")[0]] += 1
            
            mapping_rows.append({
                "variant_id": variant.get("stable_variant_id", ""),
                "gene": gene,
                "UniProt_accession": structure_data.get("uniprot_id", ""),
                "protein_position": position if position else "",
                "WT_AA": ref_aa,
                "mutant_AA": alt_aa,
                "structure_file": structure_data.get("structure_file", ""),
                "mapping_status": "success" if is_valid else "failed",
                "failure_reason": reason if not is_valid else "",
            })
    
    print(f"  Total mapped: {len(mapping_rows)}")
    print(f"  Valid: {valid_count}")
    print(f"  Invalid: {len(mapping_rows) - valid_count}")
    
    # Write mapping table
    print("\nWriting output files...")
    mapping_fields = [
        "variant_id", "gene", "UniProt_accession", "protein_position",
        "WT_AA", "mutant_AA", "structure_file", "mapping_status", "failure_reason"
    ]
    
    with mapping_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=mapping_fields)
        writer.writeheader()
        writer.writerows(mapping_rows)
    
    print(f"  Mapping: {mapping_path}")
    
    # Validation report
    validation = {
        "passed": failed == 0 and valid_count >= len(cohort_data) * 0.90,
        "checks": {
            "structures_acquired": {
                "expected": len(genes_data),
                "actual": successful,
                "passed": successful == len(genes_data),
            },
            "variant_mappings": {
                "total": len(mapping_rows),
                "valid": valid_count,
                "invalid": len(mapping_rows) - valid_count,
                "success_rate": round(valid_count / len(mapping_rows) * 100, 2) if mapping_rows else 0,
                "passed": valid_count >= len(mapping_rows) * 0.90,
            },
            "invalid_reasons": dict(invalid_reasons),
        }
    }
    
    validation_path.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    print(f"  Validation: {validation_path}")
    
    # Generate report
    report = {
        "phase": "Phase 3 - AlphaFold Structure Integration",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source": "AlphaFold Protein Structure Database",
        "source_api": ALPHAFOLD_API,
        "input": {
            "genes_file": str(genes_path),
            "cohort_file": str(cohort_path),
            "gene_count": len(genes_data),
            "variant_count": len(cohort_data),
        },
        "structures": {
            "total_genes": len(genes_data),
            "successful": successful,
            "failed": failed,
            "success_rate": round(successful / len(genes_data) * 100, 2) if genes_data else 0,
            "failures": [
                {"gene": r["gene"], "uniprot_id": r["uniprot_id"], "error": r.get("error", "")}
                for r in structure_results.values()
                if r["status"] != "success"
            ],
        },
        "mapping": {
            "total_variants": len(mapping_rows),
            "valid_mappings": valid_count,
            "invalid_mappings": len(mapping_rows) - valid_count,
            "success_rate": round(valid_count / len(mapping_rows) * 100, 2) if mapping_rows else 0,
            "invalid_reasons": dict(invalid_reasons),
        },
        "structure_details": {
            gene: {
                "uniprot_id": data.get("uniprot_id"),
                "status": data.get("status"),
                "protein_length": data.get("protein_length_cif"),
                "file": str(data.get("structure_file", "")),
            }
            for gene, data in structure_results.items()
        },
        "output_files": {
            "structures_directory": str(structure_dir),
            "mapping_file": str(mapping_path),
            "report_file": str(report_path),
            "validation_file": str(validation_path),
        },
        "storage_location": str(structure_dir),
        "storage_note": "Large structure files are stored outside GitHub in data/raw/alphafold_phase3/",
    }
    
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"  Report: {report_path}")
    
    print("\n" + "=" * 80)
    print("PHASE 3 COMPLETE")
    print("=" * 80)
    
    return report


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Phase 3: AlphaFold Structure Integration")
    parser.add_argument("--genes", type=Path, default=DEFAULT_GENES)
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--structure-dir", type=Path, default=DEFAULT_STRUCTURE_DIR)
    parser.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--validation", type=Path, default=DEFAULT_VALIDATION)
    
    args = parser.parse_args()
    
    report = integrate_alphafold(
        genes_path=args.genes,
        cohort_path=args.cohort,
        structure_dir=args.structure_dir,
        mapping_path=args.mapping,
        report_path=args.report,
        validation_path=args.validation,
    )
    
    print(f"\nStructures acquired: {report['structures']['successful']}/{report['structures']['total_genes']}")
    print(f"Valid mappings: {report['mapping']['valid_mappings']}/{report['mapping']['total_variants']} ({report['mapping']['success_rate']}%)")


if __name__ == "__main__":
    main()
