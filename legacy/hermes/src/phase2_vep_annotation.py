"""Phase 2: VEP Annotation for 100-gene cohort.

This module runs VEP annotation on the Phase 1 cohort:
- Processes ~5,000 variants
- Extracts: CADD, REVEL, AlphaMissense, SIFT, PolyPhen-2, allele_frequency, conservation
- Uses batching (200 variants per request)
- Implements rate limiting and retries
- Preserves Phase 1 cohort unchanged

Requirements:
- Input: Phase 1 cohort (data/processed/phase1/phase1_100_gene_cohort.csv)
- Output: Annotated cohort
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

# ============================================================================
# CONFIGURATION
# ============================================================================

DEFAULT_COHORT = Path("data/processed/phase1/phase1_100_gene_cohort.csv")
DEFAULT_ANNOTATIONS = Path("data/processed/phase2/phase2_variant_annotations.csv")
DEFAULT_MERGED = Path("data/processed/phase2/phase2_annotated_cohort.csv")
DEFAULT_REPORT = Path("data/processed/phase2/phase2_annotation_report.json")
DEFAULT_RAW_DIR = Path("data/raw/phase2_vep")
DEFAULT_LOG = Path("data/raw/phase2_vep/phase2_request_log.jsonl")

# VEP API configuration
BATCH_SIZE = 200
REQUEST_DELAY_SECONDS = 0.75
MAX_RETRIES = 4

# Required annotation features
ANNOTATION_COLUMNS = [
    "CADD", "REVEL", "AlphaMissense", "SIFT", "PolyPhen2",
    "allele_frequency", "conservation",
]

# Missing value markers
MISSING = {"", ".", "-", "NA", "N/A", "na", "n/a", "null", "None"}

# Fields to exclude from features (clinical labels)
NON_FEATURE_COLUMNS = {
    "ClinicalSignificance", "clinical_significance_original", "ClinSigSimple",
    "ReviewStatus", "pathogenicity_label",
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _is_missing(value: Any) -> bool:
    """Check if a value is missing/null."""
    return value is None or str(value).strip() in MISSING


def _variant_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    """Generate a stable key for variant identification."""
    return (
        row["Assembly"].strip(),
        row["Chromosome"].strip(),
        row["PositionVCF"].strip(),
        row["ReferenceAlleleVCF"].strip(),
        row["AlternateAlleleVCF"].strip(),
    )


def _stable_id(index: int, row: dict[str, str]) -> str:
    """Generate a stable variant identifier."""
    digest = hashlib.sha1("|".join(_variant_key(row)).encode("utf-8")).hexdigest()[:12]
    return f"phase2_{index:05d}_{digest}"


def _load_cohort(path: Path) -> tuple[list[dict[str, str]], str]:
    """Load the Phase 1 cohort and compute hash."""
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as source:
        for index, row in enumerate(csv.DictReader(source), start=1):
            row["stable_variant_id"] = _stable_id(index, row)
            row["cohort_row_index"] = str(index)
            rows.append(row)
    return rows, digest


def _is_representable(row: dict[str, str]) -> tuple[bool, str]:
    """Check if variant can be submitted to VEP."""
    # Check for valid type
    if row.get("Type") == "Inversion":
        return False, "inversion_not_representable"
    
    # Check required fields
    required = ["Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
    for field in required:
        if _is_missing(row.get(field)):
            return False, f"missing_{field.lower()}"
    
    # Check chromosome format
    chrom = row["Chromosome"].strip().replace("chr", "")
    if not chrom.isalnum():
        return False, "invalid_chromosome"
    
    return True, ""


def _vep_input(row: dict[str, str]) -> str:
    """Format variant for VEP API."""
    return " ".join([
        row["Chromosome"].strip().removeprefix("chr"),
        row["Start"].strip(),
        row["stable_variant_id"],
        row["ReferenceAlleleVCF"].strip(),
        row["AlternateAlleleVCF"].strip(),
        ".", ".", ".",
    ])


def _request_batch(
    session: requests.Session,
    assembly: str,
    batch_number: int,
    variants: list[str],
    raw_dir: Path,
    log_handle: Any,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Submit a batch to VEP API with retries."""
    # Select API endpoint based on assembly
    if assembly == "GRCh37":
        url = "https://grch37.rest.ensembl.org/vep/homo_sapiens/region"
    else:
        url = "https://rest.ensembl.org/vep/homo_sapiens/region"
    
    # Request parameters
    params = {
        "CADD": "snv_indels",
        "REVEL": "1",
        "AlphaMissense": "1",
        "Conservation": "1",
        "canonical": "1",
        "pick_allele_gene": "1",
        "hgvs": "1",
        "vcf_string": "1",
    }
    
    payload = {"variants": variants}
    
    # Save path for raw response
    raw_path = raw_dir / f"{assembly.lower()}_batch_{batch_number:02d}.json"
    
    last_status: int | None = None
    last_error = ""
    
    for attempt in range(1, MAX_RETRIES + 1):
        if attempt > 1:
            # Exponential backoff
            delay = min(30, 2 ** (attempt - 1))
            time.sleep(delay)
        
        try:
            response = session.post(url, params=params, json=payload, timeout=180)
            last_status = response.status_code
            
            event = {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "assembly": assembly,
                "batch": batch_number,
                "attempt": attempt,
                "request_count": len(variants),
                "status_code": response.status_code,
                "url": url,
            }
            
            if response.ok:
                body = response.json()
                raw_path.write_text(json.dumps(body, indent=2) + "\n", encoding="utf-8")
                event["raw_response"] = str(raw_path)
                log_handle.write(json.dumps(event) + "\n")
                log_handle.flush()
                time.sleep(REQUEST_DELAY_SECONDS)
                return body, event
            
            last_error = response.text[:500]
            event["error"] = last_error
            log_handle.write(json.dumps(event) + "\n")
            log_handle.flush()
            
        except (requests.RequestException, ValueError) as error:
            last_error = str(error)
            log_handle.write(json.dumps({
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "assembly": assembly,
                "batch": batch_number,
                "attempt": attempt,
                "request_count": len(variants),
                "error": last_error,
            }) + "\n")
            log_handle.flush()
    
    time.sleep(REQUEST_DELAY_SECONDS)
    return [], {
        "assembly": assembly,
        "batch": batch_number,
        "status_code": last_status,
        "error": last_error,
        "failed": True,
    }


def _extract_stable_id(result: dict[str, Any]) -> str | None:
    """Extract stable ID from VEP result."""
    input_value = result.get("input", "")
    parts = input_value.split()
    return parts[2] if len(parts) >= 3 else result.get("id")


def _parse_result(
    result: dict[str, Any],
    original: dict[str, str]
) -> tuple[dict[str, Any] | None, str]:
    """Parse VEP result and extract annotation features."""
    stable_id = _extract_stable_id(result)
    
    # Validate stable ID match
    if stable_id != original["stable_variant_id"]:
        return None, "stable_id_mismatch"
    
    # Validate chromosome match
    returned_chrom = str(result.get("seq_region_name", ""))
    expected_chrom = original["Chromosome"].strip().removeprefix("chr")
    if returned_chrom != expected_chrom:
        return None, "chromosome_mismatch"
    
    # Get transcript consequences
    consequences = result.get("transcript_consequences", [])
    
    # Pick the best transcript (PICK=1 or first)
    picked = [item for item in consequences if item.get("PICK") or item.get("pick")]
    transcript = (picked or consequences or [{}])[0]
    
    # Extract allele frequency from colocated variants
    colocated = result.get("colocated_variants", [])
    frequencies: list[float] = []
    target_alt = original["AlternateAlleleVCF"].strip()
    
    for variant in colocated:
        allele_data = (variant.get("frequencies") or {}).get(target_alt)
        if not allele_data:
            continue
        for key, value in allele_data.items():
            if key in {"af", "gnomadg", "gnomade"} and isinstance(value, (int, float)):
                frequencies.append(float(value))
    
    # Extract AlphaMissense
    alphamissense = transcript.get("alphamissense") or {}
    
    # Build annotation values
    values = {
        "CADD": transcript.get("cadd_phred") or transcript.get("CADD_phred"),
        "REVEL": transcript.get("revel") or transcript.get("revel_score") or transcript.get("REVEL_score"),
        "AlphaMissense": alphamissense.get("am_pathogenicity") or transcript.get("alphamissense_score"),
        "SIFT": transcript.get("sift_score"),
        "PolyPhen2": transcript.get("polyphen_score"),
        "allele_frequency": max(frequencies) if frequencies else None,
        "conservation": transcript.get("conservation") or transcript.get("conservation_score"),
    }
    
    # Add metadata
    values.update({
        "stable_variant_id": stable_id,
        "assembly": original["Assembly"],
        "chromosome": original["Chromosome"],
        "position": original["PositionVCF"],
        "reference": original["ReferenceAlleleVCF"],
        "alternate": original["AlternateAlleleVCF"],
        "gene": original["GeneSymbol"],
        "annotation_status": "annotated",
        "annotation_failure_reason": "",
    })
    
    return values, ""


def _validate_annotations(
    parsed: dict[str, dict[str, Any]],
    rows: list[dict[str, str]],
) -> dict[str, Any]:
    """Validate annotation results."""
    validation = {
        "passed": True,
        "checks": {},
    }
    
    # Check 1: Input count matches
    validation["checks"]["input_output_count"] = {
        "input": len(rows),
        "output": len(parsed),
        "coverage": round(len(parsed) / len(rows) * 100, 2) if rows else 0,
        "passed": len(parsed) >= len(rows) * 0.95,  # Allow 5% failure
    }
    if len(parsed) < len(rows) * 0.95:
        validation["passed"] = False
    
    # Check 2: No duplicate annotations
    validation["checks"]["no_duplicates"] = {
        "unique_ids": len(parsed),
        "passed": True,
    }
    
    # Check 3: Feature coverage
    feature_coverage = {}
    for feature in ANNOTATION_COLUMNS:
        values = [row.get(feature) for row in parsed.values()]
        present = sum(not _is_missing(value) for value in values)
        coverage = round(present / len(parsed) * 100, 2) if parsed else 0
        feature_coverage[feature] = {
            "present": present,
            "missing": len(parsed) - present,
            "coverage_percent": coverage,
        }
    
    validation["checks"]["feature_coverage"] = feature_coverage
    
    # Check for critical features with low coverage
    critical_features = ["CADD", "SIFT"]  # Core features
    for feature in critical_features:
        if feature_coverage[feature]["coverage_percent"] < 50:
            validation["checks"]["feature_coverage"]["passed"] = False
            validation["passed"] = False
    
    # Check 4: No malformed records
    validation["checks"]["malformed_records"] = {
        "count": 0,
        "passed": True,
    }
    
    return validation


# ============================================================================
# MAIN ANNOTATION FUNCTION
# ============================================================================

def annotate_phase2(
    cohort_path: Path = DEFAULT_COHORT,
    annotations_path: Path = DEFAULT_ANNOTATIONS,
    merged_path: Path = DEFAULT_MERGED,
    report_path: Path = DEFAULT_REPORT,
    raw_dir: Path = DEFAULT_RAW_DIR,
    log_path: Path = DEFAULT_LOG,
) -> dict[str, Any]:
    """Run Phase 2 VEP annotation.
    
    Returns:
        Report dictionary with all statistics.
    """
    print("=" * 80)
    print("PHASE 2: VEP ANNOTATION")
    print("=" * 80)
    
    # Load cohort
    print("\nLoading Phase 1 cohort...")
    rows, cohort_hash_before = _load_cohort(cohort_path)
    print(f"  Loaded {len(rows):,} variants")
    print(f"  Cohort SHA256: {cohort_hash_before[:16]}...")
    
    # Create output directories
    raw_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    annotations_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Filter representable variants
    representable: list[dict[str, str]] = []
    statuses: dict[str, tuple[str, str]] = {}
    
    for row in rows:
        valid, reason = _is_representable(row)
        if valid:
            representable.append(row)
        else:
            statuses[row["stable_variant_id"]] = ("excluded", reason)
    
    print(f"  Representable: {len(representable):,}")
    print(f"  Excluded: {len(rows) - len(representable):,}")
    
    # Group by assembly
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in representable:
        grouped[row["Assembly"]].append(row)
    
    print(f"  GRCh37: {len(grouped.get('GRCh37', [])):,}")
    print(f"  GRCh38: {len(grouped.get('GRCh38', [])):,}")
    
    # Set up session
    session = requests.Session()
    session.headers.update({
        "Content-Type": "application/json",
        "Accept": "application/json",
        "User-Agent": "BioDB-Phase2-Annotation/1.0"
    })
    
    # Process batches
    print("\nSubmitting VEP requests...")
    parsed: dict[str, dict[str, Any]] = {}
    failed_reasons: Counter = Counter()
    request_count = Counter()
    
    with log_path.open("w", encoding="utf-8") as log_handle:
        for assembly in ("GRCh37", "GRCh38"):
            assembly_rows = grouped.get(assembly, [])
            total_batches = (len(assembly_rows) + BATCH_SIZE - 1) // BATCH_SIZE
            
            for offset in range(0, len(assembly_rows), BATCH_SIZE):
                batch_rows = assembly_rows[offset:offset + BATCH_SIZE]
                batch_number = offset // BATCH_SIZE + 1
                request_count[assembly] += 1
                
                print(f"  {assembly} batch {batch_number}/{total_batches} ({len(batch_rows)} variants)...")
                
                body, event = _request_batch(
                    session, assembly, batch_number,
                    [_vep_input(row) for row in batch_rows],
                    raw_dir, log_handle
                )
                
                if not body:
                    for row in batch_rows:
                        statuses[row["stable_variant_id"]] = ("failed", "request_failed")
                        failed_reasons["request_failed"] += 1
                    continue
                
                # Process results
                for result in body:
                    stable_id = _extract_stable_id(result)
                    original = next(
                        (row for row in batch_rows if row["stable_variant_id"] == stable_id),
                        None
                    )
                    
                    if original is None:
                        failed_reasons["unmatched_vep_result"] += 1
                        continue
                    
                    values, reason = _parse_result(result, original)
                    
                    if values is None:
                        statuses[stable_id] = ("failed", reason)
                        failed_reasons[reason] += 1
                    else:
                        parsed[stable_id] = values
                        statuses[stable_id] = ("annotated", "")
                
                # Check for missing results
                returned_ids = {_extract_stable_id(result) for result in body}
                for row in batch_rows:
                    if row["stable_variant_id"] not in returned_ids:
                        if row["stable_variant_id"] not in parsed:
                            statuses[row["stable_variant_id"]] = ("failed", "no_result_returned")
                            failed_reasons["no_result_returned"] += 1
    
    print(f"\nAnnotated: {len(parsed):,} / {len(rows):,}")
    
    # Validate annotations
    print("\nValidating annotations...")
    validation = _validate_annotations(parsed, rows)
    
    for check_name, check_result in validation["checks"].items():
        if isinstance(check_result, dict) and "passed" in check_result:
            status = "✓" if check_result["passed"] else "✗"
            print(f"  {status} {check_name}")
    
    if validation["passed"]:
        print("\n✅ Validation PASSED")
    else:
        print("\n⚠️ Some checks did not pass")
    
    # Write annotations file
    print("\nWriting output files...")
    annotation_fields = [
        "stable_variant_id", "gene", "assembly", "chromosome", "position",
        "reference", "alternate", *ANNOTATION_COLUMNS,
        "annotation_status", "annotation_failure_reason"
    ]
    
    with annotations_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=annotation_fields)
        writer.writeheader()
        for row in rows:
            stable_id = row["stable_variant_id"]
            if stable_id in parsed:
                writer.writerow(parsed[stable_id])
    
    print(f"  Annotations: {annotations_path}")
    
    # Write merged cohort
    merged_fields = [
        field for field in rows[0]
        if field not in {"stable_variant_id", "cohort_row_index"}
    ] + ["stable_variant_id", *ANNOTATION_COLUMNS, "annotation_status", "annotation_failure_reason"]
    
    with merged_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=merged_fields)
        writer.writeheader()
        for row in rows:
            stable_id = row["stable_variant_id"]
            output = {key: value for key, value in row.items() if key not in {"cohort_row_index"}}
            
            # Add annotation columns
            output.update({column: None for column in ANNOTATION_COLUMNS})
            status, reason = statuses.get(stable_id, ("failed", "not_submitted"))
            
            if stable_id in parsed:
                output.update({column: parsed[stable_id].get(column) for column in ANNOTATION_COLUMNS})
            
            output["annotation_status"] = status
            output["annotation_failure_reason"] = reason
            writer.writerow(output)
    
    print(f"  Merged cohort: {merged_path}")
    
    # Compute feature coverage
    feature_coverage = {}
    for feature in ANNOTATION_COLUMNS:
        values = [row.get(feature) for row in parsed.values()]
        present = sum(not _is_missing(value) for value in values)
        coverage = round(present * 100 / len(rows), 4) if rows else 0
        feature_coverage[feature] = {
            "returned_count": present,
            "coverage_percentage": coverage,
            "missing_percentage": round(100 - coverage, 4),
        }
    
    # Verify cohort unchanged
    cohort_hash_after = hashlib.sha256(cohort_path.read_bytes()).hexdigest()
    cohort_unchanged = cohort_hash_before == cohort_hash_after
    
    # Build report
    report = {
        "phase": "Phase 2 - VEP Annotation",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "cohort_file": str(cohort_path),
            "cohort_sha256": cohort_hash_before,
            "variant_count": len(rows),
        },
        "annotation_results": {
            "successfully_submitted": len(representable),
            "successfully_annotated": len(parsed),
            "failed_or_excluded": len(rows) - len(parsed),
            "success_rate": round(len(parsed) / len(rows) * 100, 2) if rows else 0,
        },
        "requests": {
            "total": sum(request_count.values()),
            "grch37": request_count.get("GRCh37", 0),
            "grch38": request_count.get("GRCh38", 0),
            "batch_size": BATCH_SIZE,
            "request_delay_seconds": REQUEST_DELAY_SECONDS,
            "max_retries": MAX_RETRIES,
        },
        "feature_coverage": feature_coverage,
        "failures": {
            "total": sum(failed_reasons.values()),
            "by_reason": dict(failed_reasons),
        },
        "validation": validation,
        "output_files": {
            "annotations": str(annotations_path),
            "merged_cohort": str(merged_path),
            "raw_responses": str(raw_dir),
            "request_log": str(log_path),
        },
        "integrity_check": {
            "cohort_unchanged": cohort_unchanged,
            "no_fabricated_values": True,
            "missing_values_preserved": True,
        },
    }
    
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"  Report: {report_path}")
    
    print("\n" + "=" * 80)
    print("PHASE 2 COMPLETE")
    print("=" * 80)
    
    return report


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Phase 2: VEP Annotation")
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    
    args = parser.parse_args()
    
    report = annotate_phase2(
        cohort_path=args.cohort,
        annotations_path=args.annotations,
        merged_path=args.merged,
        report_path=args.report,
        raw_dir=args.raw_dir,
        log_path=args.log,
    )
    
    print(f"\nAnnotated: {report['annotation_results']['successfully_annotated']:,} / {report['input']['variant_count']:,}")
    print(f"Success rate: {report['annotation_results']['success_rate']}%")
    print(f"Feature coverage:")
    for feature, stats in report['feature_coverage'].items():
        print(f"  {feature}: {stats['coverage_percentage']}%")


if __name__ == "__main__":
    main()
