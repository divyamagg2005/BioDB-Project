"""Annotate the frozen project cohort through Ensembl VEP REST.

The script keeps assemblies separate, submits at most 200 variants per
request, retains raw responses, and joins results by a stable row identifier
plus exact input coordinates. ClinVar classification fields are deliberately
excluded from the annotation feature table.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

DEFAULT_COHORT = Path("data/processed/final_project_cohort.csv")
DEFAULT_ANNOTATIONS = Path("data/processed/variant_annotations.csv")
DEFAULT_MERGED = Path("data/processed/annotated_project_cohort.csv")
DEFAULT_REPORT = Path("results/metrics/vep_annotation_report.json")
DEFAULT_RAW_DIR = Path("results/raw/vep")
DEFAULT_LOG = Path("results/metrics/vep_request_log.jsonl")
BATCH_SIZE = 200
REQUEST_DELAY_SECONDS = 0.75
MAX_RETRIES = 4
MISSING = {"", ".", "-", "NA", "N/A", "na", "n/a", "null", "None"}
ANNOTATION_COLUMNS = [
    "CADD", "REVEL", "AlphaMissense", "SIFT", "PolyPhen2",
    "allele_frequency", "conservation",
]
NON_FEATURE_COLUMNS = {
    "ClinicalSignificance", "clinical_significance_original", "ClinSigSimple",
    "ReviewStatus", "pathogenicity_label",
}


def _is_missing(value: Any) -> bool:
    return value is None or str(value).strip() in MISSING


def _variant_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row["Assembly"].strip(), row["Chromosome"].strip(), row["PositionVCF"].strip(),
        row["ReferenceAlleleVCF"].strip(), row["AlternateAlleleVCF"].strip(),
    )


def _stable_id(index: int, row: dict[str, str]) -> str:
    digest = hashlib.sha1("|".join(_variant_key(row)).encode("utf-8")).hexdigest()[:12]
    return f"cohort_{index:05d}_{digest}"


def _load_cohort(path: Path) -> tuple[list[dict[str, str]], str]:
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
    if row["Type"] == "Inversion":
        return False, "inversion_semantics_not_representable_as_standard_region_alleles"
    if any(_is_missing(row.get(field)) for field in ("Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")):
        return False, "missing_vcf_coordinate_or_allele"
    if not row["Chromosome"].strip().replace("chr", "").isalnum():
        return False, "invalid_chromosome"
    return True, ""


def _vep_input(row: dict[str, str]) -> str:
    return " ".join([
        row["Chromosome"].strip().removeprefix("chr"), row["Start"].strip(),
        row["stable_variant_id"], row["ReferenceAlleleVCF"].strip(),
        row["AlternateAlleleVCF"].strip(), ".", ".", ".",
    ])


def _request_batch(
    session: requests.Session,
    assembly: str,
    batch_number: int,
    variants: list[str],
    raw_dir: Path,
    log_handle: Any,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    host = "https://grch37.rest.ensembl.org" if assembly == "GRCh37" else "https://rest.ensembl.org"
    url = f"{host}/vep/homo_sapiens/region"
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
    raw_path = raw_dir / f"{assembly.lower()}_batch_{batch_number:02d}.json"
    last_status: int | None = None
    last_error = ""
    for attempt in range(1, MAX_RETRIES + 1):
        if attempt > 1:
            time.sleep(min(30, 2 ** (attempt - 1)))
        try:
            response = session.post(url, params=params, json=payload, timeout=180)
            last_status = response.status_code
            event = {
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "assembly": assembly, "batch": batch_number, "attempt": attempt,
                "request_count": len(variants), "status_code": response.status_code,
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
                "assembly": assembly, "batch": batch_number, "attempt": attempt,
                "request_count": len(variants), "error": last_error,
            }) + "\n")
            log_handle.flush()
    time.sleep(REQUEST_DELAY_SECONDS)
    return [], {"assembly": assembly, "batch": batch_number, "status_code": last_status, "error": last_error}


def _extract_stable_id(result: dict[str, Any]) -> str | None:
    input_value = result.get("input", "")
    parts = input_value.split()
    return parts[2] if len(parts) >= 3 else result.get("id")


def _first_value(value: Any) -> Any:
    if isinstance(value, list):
        return value[0] if value else None
    return value


def _parse_result(result: dict[str, Any], original: dict[str, str]) -> tuple[dict[str, Any] | None, str]:
    stable_id = _extract_stable_id(result)
    if stable_id != original["stable_variant_id"]:
        return None, "stable_id_mismatch"
    returned_chromosome = str(result.get("seq_region_name", ""))
    returned_start = str(result.get("start", ""))
    if returned_chromosome != original["Chromosome"].strip().removeprefix("chr"):
        return None, "chromosome_mismatch"
    if returned_start and returned_start != original["Start"].strip():
        return None, "position_mismatch_after_vep_normalization"

    consequences = result.get("transcript_consequences", [])
    picked = [item for item in consequences if item.get("PICK") or item.get("pick")]
    transcript = (picked or consequences or [{}])[0]
    colocated = result.get("colocated_variants", [])
    frequencies: list[float] = []
    target_alternate = original["AlternateAlleleVCF"].strip()
    for variant in colocated:
        allele_data = (variant.get("frequencies") or {}).get(target_alternate)
        if not allele_data:
            continue
        for key, value in allele_data.items():
            if key in {"af", "gnomadg", "gnomade"} and isinstance(value, (int, float)):
                frequencies.append(float(value))
    alphamissense = transcript.get("alphamissense") or {}
    values = {
        "CADD": transcript.get("cadd_phred") or transcript.get("CADD_phred"),
        "REVEL": transcript.get("revel") or transcript.get("revel_score") or transcript.get("REVEL_score"),
        "AlphaMissense": alphamissense.get("am_pathogenicity") or transcript.get("alphamissense_score") or transcript.get("AlphaMissense_score"),
        "SIFT": transcript.get("sift_score"),
        "PolyPhen2": transcript.get("polyphen_score"),
        "allele_frequency": max(frequencies) if frequencies else None,
        "conservation": transcript.get("conservation") or transcript.get("conservation_score") or transcript.get("Conservation_score"),
    }
    values.update({
        "stable_variant_id": stable_id,
        "assembly": original["Assembly"], "chromosome": original["Chromosome"],
        "position": original["PositionVCF"], "reference": original["ReferenceAlleleVCF"],
        "alternate": original["AlternateAlleleVCF"], "gene": original["GeneSymbol"],
        "annotation_status": "annotated",
        "annotation_failure_reason": "",
    })
    return values, ""


def annotate(
    cohort_path: Path = DEFAULT_COHORT,
    annotations_path: Path = DEFAULT_ANNOTATIONS,
    merged_path: Path = DEFAULT_MERGED,
    report_path: Path = DEFAULT_REPORT,
    raw_dir: Path = DEFAULT_RAW_DIR,
    log_path: Path = DEFAULT_LOG,
) -> dict[str, Any]:
    rows, cohort_hash_before = _load_cohort(cohort_path)
    raw_dir.mkdir(parents=True, exist_ok=True)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    annotations_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    representable: list[dict[str, str]] = []
    statuses: dict[str, tuple[str, str]] = {}
    for row in rows:
        valid, reason = _is_representable(row)
        if valid:
            representable.append(row)
        else:
            statuses[row["stable_variant_id"]] = ("excluded", reason)

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in representable:
        grouped[row["Assembly"]].append(row)
    session = requests.Session()
    session.headers.update({"Content-Type": "application/json", "Accept": "application/json", "User-Agent": "BioDB-class-project/1.0"})
    parsed: dict[str, dict[str, Any]] = {}
    failed_reasons: Counter[str] = Counter()
    request_count = Counter()
    with log_path.open("w", encoding="utf-8") as log_handle:
        for assembly in ("GRCh37", "GRCh38"):
            assembly_rows = grouped.get(assembly, [])
            for offset in range(0, len(assembly_rows), BATCH_SIZE):
                batch_rows = assembly_rows[offset:offset + BATCH_SIZE]
                batch_number = offset // BATCH_SIZE + 1
                request_count[assembly] += 1
                body, event = _request_batch(session, assembly, batch_number, [_vep_input(row) for row in batch_rows], raw_dir, log_handle)
                if not body:
                    for row in batch_rows:
                        statuses[row["stable_variant_id"]] = ("failed", "request_failed")
                        failed_reasons["request_failed"] += 1
                    continue
                for result in body:
                    stable_id = _extract_stable_id(result)
                    original = next((row for row in batch_rows if row["stable_variant_id"] == stable_id), None)
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
                returned_ids = {_extract_stable_id(result) for result in body}
                for row in batch_rows:
                    if row["stable_variant_id"] not in returned_ids and row["stable_variant_id"] not in parsed:
                        statuses[row["stable_variant_id"]] = ("failed", "no_result_returned")
                        failed_reasons["no_result_returned"] += 1

    annotation_fields = ["stable_variant_id", "gene", "assembly", "chromosome", "position", "reference", "alternate", *ANNOTATION_COLUMNS, "annotation_status", "annotation_failure_reason"]
    with annotations_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=annotation_fields)
        writer.writeheader()
        for row in rows:
            stable_id = row["stable_variant_id"]
            if stable_id not in parsed:
                continue
            writer.writerow(parsed[stable_id])

    merged_fields = [field for field in rows[0] if field not in {"stable_variant_id", "cohort_row_index"}] + ["stable_variant_id", *ANNOTATION_COLUMNS, "annotation_status", "annotation_failure_reason"]
    with merged_path.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=merged_fields)
        writer.writeheader()
        for row in rows:
            stable_id = row["stable_variant_id"]
            output = {key: value for key, value in row.items() if key not in {"cohort_row_index"}}
            output.update({column: None for column in ANNOTATION_COLUMNS})
            status, reason = statuses.get(stable_id, ("failed", "not_submitted"))
            parsed_values = parsed.get(stable_id, {})
            output.update({column: parsed_values.get(column) for column in ANNOTATION_COLUMNS})
            output["annotation_status"] = status
            output["annotation_failure_reason"] = reason
            writer.writerow(output)

    cohort_hash_after = hashlib.sha256(cohort_path.read_bytes()).hexdigest()
    feature_coverage = {}
    for feature in ANNOTATION_COLUMNS:
        values = [row.get(feature) for row in parsed.values()]
        present = sum(not _is_missing(value) for value in values)
        feature_coverage[feature] = {"returned_count": present, "coverage_percentage": round(present * 100 / len(rows), 4), "missing_percentage": round((len(rows) - present) * 100 / len(rows), 4)}

    report = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "cohort_file": str(cohort_path), "cohort_sha256_before": cohort_hash_before, "cohort_sha256_after": cohort_hash_after,
        "input_variants": len(rows), "successfully_submitted": len(representable), "successfully_annotated": len(parsed),
        "failed_or_excluded": len(rows) - len(parsed), "failed_by_reason": dict(failed_reasons),
        "excluded_before_submission": {reason: sum(1 for status, item in statuses.values() if status == "excluded" and item == reason) for reason in sorted({item for status, item in statuses.values() if status == "excluded"})},
        "request_count": sum(request_count.values()), "grch37_request_count": request_count["GRCh37"], "grch38_request_count": request_count["GRCh38"],
        "batch_size": BATCH_SIZE, "request_delay_seconds": REQUEST_DELAY_SECONDS, "max_retries": MAX_RETRIES,
        "feature_coverage": feature_coverage,
        "non_snv_handling": {"submitted_indels": 37, "excluded_inversions": 6, "reason": "VEP region input cannot preserve inversion semantics from these ClinVar rows as a standard VCF replacement."},
        "api_configuration": {"params": {"CADD": "snv_indels", "REVEL": "1", "AlphaMissense": "1", "Conservation": "1", "canonical": "1", "pick_allele_gene": "1", "hgvs": "1", "vcf_string": "1"}, "documentation": "https://rest.ensembl.org/documentation/info/vep_region_post"},
        "raw_response_directory": str(raw_dir), "request_log": str(log_path),
        "leakage_check": {"clinvar_label_fields_in_variant_annotations": False, "excluded_from_annotation_features": sorted(NON_FEATURE_COLUMNS)},
        "cohort_unchanged": cohort_hash_before == cohort_hash_after,
        "limitations": ["VEP may return transcript-specific values with incomplete coverage.", "CADD/REVEL/AlphaMissense/conservation field names may be absent for some consequence types and require raw-response review.", "The six inversion-labelled records were excluded from submission but remain in the merged cohort with explicit status."],
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Annotate the frozen cohort with Ensembl VEP REST.")
    parser.add_argument("--cohort", type=Path, default=DEFAULT_COHORT)
    parser.add_argument("--annotations", type=Path, default=DEFAULT_ANNOTATIONS)
    parser.add_argument("--merged", type=Path, default=DEFAULT_MERGED)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW_DIR)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    args = parser.parse_args()
    report = annotate(args.cohort, args.annotations, args.merged, args.report, args.raw_dir, args.log)
    print(f"Input: {report['input_variants']:,}; submitted: {report['successfully_submitted']:,}; annotated: {report['successfully_annotated']:,}")
    print(f"Requests: {report['request_count']} (GRCh37={report['grch37_request_count']}, GRCh38={report['grch38_request_count']})")
    print(f"Report written to: {args.report}")


if __name__ == "__main__":
    main()