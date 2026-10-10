"""Phase 2: annotate the Phase 1 cohort with Ensembl VEP REST.

Each cohort variant is submitted as a GRCh38 VCF-style record (``PositionVCF``
with the VCF alleles) in batches of 200. Raw responses are saved per batch in
``data/raw/phase2_vep/`` and reused on reruns, so an interrupted run resumes
where it stopped and a finished run can be re-parsed offline.

From each response one transcript is chosen: a transcript of the cohort gene
whose protein position and amino-acid change equal the ClinVar protein change,
preferring MANE Select, then the Ensembl canonical transcript. Features are
taken from that transcript, plus gnomAD allele frequencies for the exact
alternate allele.

ClinVar fields that VEP returns (``clin_sig``, ``var_synonyms``, phenotype
flags) are never read: they would leak the label.

Run with ``python -m biodb.phase2_vep``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

from biodb import config
from biodb.common import sha256_file, write_json

VEP_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
VEP_PARAMS = {
    "CADD": "1",
    "REVEL": "1",
    "AlphaMissense": "1",
    "Conservation": "1",
    "canonical": "1",
    "mane": "1",
    "protein": "1",
    "uniprot": "1",
    "numbers": "1",
    "variant_class": "1",
    "af_gnomade": "1",
    "af_gnomadg": "1",
}
BATCH_SIZE = 200
REQUEST_DELAY = 1.0
MAX_ATTEMPTS = 5

FEATURES = [
    "cadd_phred", "cadd_raw", "revel", "alphamissense", "sift_score", "polyphen_score",
    "gerp_conservation", "gnomade_af", "gnomadg_af", "gnomad_af_max",
]
OUTPUT_FIELDS = [
    "variant_id", "annotation_status", "failure_reason",
    "vep_transcript", "vep_protein_id", "transcript_choice", "mane_select", "uniprot_isoform",
    "consequence", "impact", "vep_protein_position", "vep_amino_acids",
    *FEATURES, "scores_from_other_transcript",
    "alphamissense_class", "sift_prediction", "polyphen_prediction",
]


# ---------------------------------------------------------------------------
# Request handling
# ---------------------------------------------------------------------------

def vep_input(row: dict[str, str]) -> str:
    """VCF-style VEP region input; the ID column carries the cohort variant_id."""
    return " ".join([
        row["Chromosome"], row["PositionVCF"], row["variant_id"],
        row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"], ".", ".", ".",
    ])


def _chromosome_order(chromosome: str) -> tuple[int, str]:
    return (int(chromosome), "") if chromosome.isdigit() else (100, chromosome)


def batch_lines(batch: list[dict[str, str]]) -> list[str]:
    """VEP input lines for a batch, sorted by chromosome and position as VEP requires."""
    ordered = sorted(batch, key=lambda r: (_chromosome_order(r["Chromosome"]), int(r["PositionVCF"])))
    return [vep_input(r) for r in ordered]


def fetch_batch(session: requests.Session, lines: list[str], raw_path: Path, log: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Return the VEP response for one batch, from disk if already fetched."""
    if raw_path.exists():
        return json.loads(raw_path.read_text(encoding="utf-8"))
    last_error = ""
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            response = session.post(
                VEP_URL,
                params=VEP_PARAMS,
                json={"variants": lines},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=300,
            )
            log.append({
                "batch": raw_path.name, "attempt": attempt, "status": response.status_code,
                "time_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            })
            if response.status_code == 429:
                wait = float(response.headers.get("Retry-After", 2 ** attempt))
                time.sleep(wait)
                continue
            if not response.ok:
                log[-1]["error"] = response.text[:300]
            response.raise_for_status()
            body = response.json()
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            raw_path.write_text(json.dumps(body) + "\n", encoding="utf-8")
            time.sleep(REQUEST_DELAY)
            return body
        except (requests.RequestException, ValueError) as error:
            last_error = str(error)[:300]
            if isinstance(error, requests.HTTPError) and error.response is not None and error.response.status_code == 400:
                break
            log.append({"batch": raw_path.name, "attempt": attempt, "error": last_error})
            time.sleep(min(60, 2 ** attempt))
    raise RuntimeError(f"VEP batch {raw_path.name} failed after {MAX_ATTEMPTS} attempts: {last_error}")


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _number(value: Any) -> float | None:
    if isinstance(value, bool) or value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def choose_transcript(consequences: list[dict[str, Any]], row: dict[str, str]) -> tuple[dict[str, Any] | None, str]:
    """Pick the transcript whose protein change equals the ClinVar protein change.

    Returns (transcript, choice) where choice is 'mane_select', 'canonical' or
    'other_matching', or (None, reason) when no transcript matches.
    """
    gene_hits = [t for t in consequences if t.get("gene_symbol") == row["gene"]]
    if not gene_hits:
        return None, "no_transcript_for_gene"
    expected_aa = f"{row['ref_aa']}/{row['alt_aa']}"
    position = int(row["protein_position"])
    matching = [
        t for t in gene_hits
        if t.get("protein_start") == position
        and t.get("amino_acids") == expected_aa
        and "missense_variant" in t.get("consequence_terms", [])
    ]
    if not matching:
        return None, "no_transcript_matches_protein_change"
    for choice, test in (
        ("mane_select", lambda t: bool(t.get("mane_select"))),
        ("canonical", lambda t: t.get("canonical") == 1),
    ):
        hits = sorted((t for t in matching if test(t)), key=lambda t: t.get("transcript_id", ""))
        if hits:
            return hits[0], choice
    return sorted(matching, key=lambda t: t.get("transcript_id", ""))[0], "other_matching"


def allele_frequencies(result: dict[str, Any], alt: str) -> dict[str, float | None]:
    """gnomAD exome and genome frequencies of the exact alternate allele."""
    exome: list[float] = []
    genome: list[float] = []
    for colocated in result.get("colocated_variants", []):
        allele = (colocated.get("frequencies") or {}).get(alt) or {}
        if (value := _number(allele.get("gnomade"))) is not None:
            exome.append(value)
        if (value := _number(allele.get("gnomadg"))) is not None:
            genome.append(value)
    gnomade = max(exome) if exome else None
    gnomadg = max(genome) if genome else None
    present = [v for v in (gnomade, gnomadg) if v is not None]
    return {"gnomade_af": gnomade, "gnomadg_af": gnomadg, "gnomad_af_max": max(present) if present else None}


# Transcript-level scores that may be missing on the chosen transcript only
# because the score's source predates it. For these the median over other
# transcripts carrying the same amino-acid change is used instead.
FALLBACK_SCORES = {
    "revel": lambda t: t.get("revel"),
    "polyphen_score": lambda t: t.get("polyphen_score"),
    "sift_score": lambda t: t.get("sift_score"),
    "alphamissense": lambda t: (t.get("alphamissense") or {}).get("am_pathogenicity"),
}


def fill_missing_scores(values: dict[str, Any], consequences: list[dict[str, Any]], row: dict[str, str]) -> list[str]:
    """Fill FALLBACK_SCORES absent on the chosen transcript; return the names filled."""
    expected_aa = f"{row['ref_aa']}/{row['alt_aa']}"
    peers = [t for t in consequences if t.get("gene_symbol") == row["gene"] and t.get("amino_acids") == expected_aa]
    filled = []
    for name, getter in FALLBACK_SCORES.items():
        if values[name] is not None:
            continue
        found = sorted(v for t in peers if (v := _number(getter(t))) is not None)
        if found:
            middle = len(found) // 2
            values[name] = found[middle] if len(found) % 2 else (found[middle - 1] + found[middle]) / 2
            filled.append(name)
    return filled


def parse_result(result: dict[str, Any] | None, row: dict[str, str]) -> dict[str, Any]:
    base: dict[str, Any] = {field: None for field in OUTPUT_FIELDS}
    base["variant_id"] = row["variant_id"]
    if result is None:
        return {**base, "annotation_status": "failed", "failure_reason": "no_vep_result"}
    if str(result.get("seq_region_name")) != row["Chromosome"] or int(result.get("start", -1)) != int(row["PositionVCF"]):
        return {**base, "annotation_status": "failed", "failure_reason": "coordinate_mismatch"}
    if result.get("allele_string") != f"{row['ReferenceAlleleVCF']}/{row['AlternateAlleleVCF']}":
        return {**base, "annotation_status": "failed", "failure_reason": "allele_mismatch"}

    transcript, choice = choose_transcript(result.get("transcript_consequences", []), row)
    if transcript is None:
        return {**base, "annotation_status": "failed", "failure_reason": choice}

    alphamissense = transcript.get("alphamissense") or {}
    values = {
        **base,
        "annotation_status": "annotated",
        "failure_reason": "",
        "vep_transcript": transcript.get("transcript_id"),
        "vep_protein_id": transcript.get("protein_id"),
        "transcript_choice": choice,
        "mane_select": transcript.get("mane_select") or "",
        "uniprot_isoform": ";".join(transcript.get("uniprot_isoform") or []),
        "consequence": "&".join(transcript.get("consequence_terms", [])),
        "impact": transcript.get("impact"),
        "vep_protein_position": transcript.get("protein_start"),
        "vep_amino_acids": transcript.get("amino_acids"),
        "cadd_phred": _number(transcript.get("cadd_phred")),
        "cadd_raw": _number(transcript.get("cadd_raw")),
        "revel": _number(transcript.get("revel")),
        "alphamissense": _number(alphamissense.get("am_pathogenicity")),
        "alphamissense_class": alphamissense.get("am_class"),
        "sift_score": _number(transcript.get("sift_score")),
        "sift_prediction": transcript.get("sift_prediction"),
        "polyphen_score": _number(transcript.get("polyphen_score")),
        "polyphen_prediction": transcript.get("polyphen_prediction"),
        "gerp_conservation": _number(transcript.get("conservation")),
        **allele_frequencies(result, row["AlternateAlleleVCF"]),
    }
    filled = fill_missing_scores(values, result.get("transcript_consequences", []), row)
    values["scores_from_other_transcript"] = ";".join(filled)
    return values


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate(annotations: list[dict[str, Any]], cohort: list[dict[str, str]], max_failure_rate: float) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: Any = "") -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    ids = [a["variant_id"] for a in annotations]
    failed = [a for a in annotations if a["annotation_status"] != "annotated"]
    check("one_row_per_cohort_variant", sorted(ids) == sorted(r["variant_id"] for r in cohort), len(ids))
    check("no_duplicate_annotations", len(ids) == len(set(ids)))
    check("failure_rate_within_limit", len(failed) / len(annotations) <= max_failure_rate,
          {"failed": len(failed), "limit": max_failure_rate})
    check("failures_have_reason", all(a["failure_reason"] for a in failed))
    annotated = [a for a in annotations if a["annotation_status"] == "annotated"]
    cohort_by_id = {r["variant_id"]: r for r in cohort}
    check("protein_change_matches_cohort", all(
        a["vep_protein_position"] == int(cohort_by_id[a["variant_id"]]["protein_position"])
        and a["vep_amino_acids"] == f"{cohort_by_id[a['variant_id']]['ref_aa']}/{cohort_by_id[a['variant_id']]['alt_aa']}"
        for a in annotated
    ))
    ranges = {
        "revel": (0, 1), "alphamissense": (0, 1), "sift_score": (0, 1), "polyphen_score": (0, 1),
        "gnomade_af": (0, 1), "gnomadg_af": (0, 1), "gnomad_af_max": (0, 1), "cadd_phred": (0, 100),
    }
    out_of_range = {
        f: sum(1 for a in annotated if a[f] is not None and not lo <= a[f] <= hi)
        for f, (lo, hi) in ranges.items()
    }
    check("feature_values_in_range", not any(out_of_range.values()), out_of_range)
    check("no_label_columns", not {"label", "ClinicalSignificance", "ReviewStatus"} & set(OUTPUT_FIELDS))
    return checks


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def annotate(
    cohort_path: Path = config.PHASE1_DIR / "cohort.csv",
    output_dir: Path = config.PHASE2_DIR,
    raw_dir: Path = config.PHASE2_RAW_DIR,
    offline: bool = False,
    max_failure_rate: float = 0.02,
) -> dict[str, Any]:
    with cohort_path.open(encoding="utf-8", newline="") as handle:
        cohort = list(csv.DictReader(handle))
    batches = [cohort[i:i + BATCH_SIZE] for i in range(0, len(cohort), BATCH_SIZE)]
    print(f"{len(cohort):,} variants in {len(batches)} batches", flush=True)

    session = requests.Session()
    log: list[dict[str, Any]] = []
    annotations: list[dict[str, Any]] = []
    for number, batch in enumerate(batches, start=1):
        raw_path = raw_dir / f"batch_{number:03d}.json"
        if offline and not raw_path.exists():
            raise RuntimeError(f"offline mode and no saved response: {raw_path}")
        try:
            body = fetch_batch(session, batch_lines(batch), raw_path, log)
        finally:
            if log:
                with (raw_dir / "request_log.jsonl").open("a", encoding="utf-8") as handle:
                    handle.writelines(json.dumps(event) + "\n" for event in log)
                log.clear()
        by_id: dict[str, dict[str, Any]] = {}
        for result in body:
            by_id.setdefault(result.get("id", ""), result)
        annotations.extend(parse_result(by_id.get(row["variant_id"]), row) for row in batch)
        print(f"  batch {number}/{len(batches)} done", flush=True)

    checks = validate(annotations, cohort, max_failure_rate)
    passed = all(c["passed"] for c in checks)

    output_dir.mkdir(parents=True, exist_ok=True)
    annotations_path = output_dir / "vep_annotations.csv"
    with annotations_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in annotations:
            writer.writerow({k: ("" if v is None else v) for k, v in row.items()})

    labels = {r["variant_id"]: r["label"] for r in cohort}
    annotated = [a for a in annotations if a["annotation_status"] == "annotated"]
    coverage = {}
    for feature in FEATURES:
        present = [a for a in annotations if a[feature] is not None]
        coverage[feature] = {
            "present": len(present),
            "coverage": round(len(present) / len(annotations), 4),
            "coverage_pathogenic": round(sum(labels[a["variant_id"]] == "1" for a in present) / sum(v == "1" for v in labels.values()), 4),
            "coverage_benign": round(sum(labels[a["variant_id"]] == "0" for a in present) / sum(v == "0" for v in labels.values()), 4),
        }
    report = {
        "phase": 2,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": {"vep_url": VEP_URL, "parameters": VEP_PARAMS, "batch_size": BATCH_SIZE},
        "input": {"cohort": cohort_path.relative_to(config.PROJECT_ROOT).as_posix(), "cohort_sha256": sha256_file(cohort_path)},
        "variants": len(annotations),
        "annotated": len(annotated),
        "failed": len(annotations) - len(annotated),
        "failure_reasons": dict(Counter(a["failure_reason"] for a in annotations if a["failure_reason"])),
        "transcript_choice": dict(Counter(a["transcript_choice"] for a in annotated)),
        "feature_coverage": coverage,
        "raw_response_sha256": {p.name: sha256_file(p) for p in sorted(raw_dir.glob("batch_*.json"))},
        "annotations_sha256": sha256_file(annotations_path),
    }
    write_json(output_dir / "annotation_report.json", report)
    write_json(output_dir / "validation_report.json", {"passed": passed, "checks": checks})

    print(f"\nAnnotated {len(annotated):,}/{len(annotations):,}")
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['check']}  {c['detail'] if not c['passed'] else ''}")
    return {"passed": passed, "report": report, "checks": checks}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--cohort", type=Path, default=config.PHASE1_DIR / "cohort.csv")
    parser.add_argument("--output-dir", type=Path, default=config.PHASE2_DIR)
    parser.add_argument("--raw-dir", type=Path, default=config.PHASE2_RAW_DIR)
    parser.add_argument("--offline", action="store_true", help="parse saved responses only")
    args = parser.parse_args(argv)
    result = annotate(args.cohort, args.output_dir, args.raw_dir, args.offline)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
