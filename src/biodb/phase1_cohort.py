"""Phase 1: build the 100-gene ClinVar missense cohort.

Pipeline:

1. Stream ClinVar ``variant_summary.txt.gz`` and keep GRCh38 germline
   single-nucleotide missense variants with a pathogenic or benign
   classification and at least one review star.
2. Deduplicate by VariationID, then by (gene, protein change). A protein change
   seen with both labels is dropped entirely.
3. For each gene with enough variants in both classes, map the gene symbol to
   one reviewed human UniProt entry and its full-length AlphaFold DB model,
   and require the AlphaFold sequence to equal the UniProt sequence.
4. Keep only variants whose reference residue matches that sequence at the
   stated position (this removes variants named on a non-canonical isoform).
5. Rank the genes that still have enough variants in both classes and take
   the top ``TARGET_GENES``. Sample ``VARIANTS_PER_CLASS`` variants per class
   per gene at random with a fixed seed.
6. Validate the cohort and write the outputs to ``data/processed/phase1/``.

Run with ``python -m biodb.phase1_cohort``. API responses are cached under
``data/raw/phase1_api_cache/`` so reruns are offline and reproducible.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

import requests

from biodb import config
from biodb.common import parse_protein_change, parse_transcript, sha256_file, usable, write_json

# ClinVar germline classifications kept, mapped to the binary target.
LABELS = {
    "Pathogenic": 1,
    "Likely pathogenic": 1,
    "Pathogenic/Likely pathogenic": 1,
    "Benign": 0,
    "Likely benign": 0,
    "Benign/Likely benign": 0,
}

# ClinVar review status to star rating. Statuses not listed here (no assertion
# criteria, no classification, conflicting) have zero stars and are excluded.
REVIEW_STARS = {
    "criteria provided, single submitter": 1,
    "criteria provided, multiple submitters, no conflicts": 2,
    "reviewed by expert panel": 3,
    "practice guideline": 4,
}
MIN_REVIEW_STARS = 1

CHROMOSOMES = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

COHORT_FIELDS = [
    "variant_id", "VariationID", "AlleleID", "gene", "GeneID", "HGNC_ID",
    "uniprot_accession", "transcript", "protein_change", "ref_aa", "protein_position", "alt_aa",
    "label", "ClinicalSignificance", "ReviewStatus", "review_stars", "NumberSubmitters",
    "LastEvaluated", "Assembly", "Chromosome", "Start", "Stop", "PositionVCF",
    "ReferenceAlleleVCF", "AlternateAlleleVCF", "rsID", "Name",
]


# ---------------------------------------------------------------------------
# ClinVar filtering
# ---------------------------------------------------------------------------

def _name_gene(name: str) -> str:
    """Gene symbol in parentheses after the transcript, e.g. 'TP53' in 'NM_000546.6(TP53):c...'."""
    start = name.find("(")
    end = name.find(")", start + 1)
    return name[start + 1:end] if start != -1 and end != -1 else ""


def filter_row(row: dict[str, str]) -> tuple[dict[str, Any] | None, str]:
    """Return (record, "") for a kept ClinVar row, or (None, reason) for a dropped one."""
    if row.get("Assembly") != config.ASSEMBLY:
        return None, "other_assembly"
    significance = (row.get("ClinicalSignificance") or "").strip()
    if significance not in LABELS:
        return None, "label_not_pathogenic_or_benign"
    if row.get("Type") != "single nucleotide variant":
        return None, "not_snv"
    if "germline" not in (row.get("OriginSimple") or ""):
        return None, "not_germline"
    review = (row.get("ReviewStatus") or "").strip()
    stars = REVIEW_STARS.get(review, 0)
    if stars < MIN_REVIEW_STARS:
        return None, "review_below_one_star"
    name = row.get("Name") or ""
    change = parse_protein_change(name)
    if change is None:
        return None, "not_missense"
    if change[0] == "M" and change[1] == 1:
        # Initiator-methionine changes abolish translation start; they are not
        # substitutions in the folded protein.
        return None, "start_loss"
    gene = (row.get("GeneSymbol") or "").strip()
    if not usable(gene) or ";" in gene or _name_gene(name) != gene:
        return None, "ambiguous_gene"
    chromosome = (row.get("Chromosome") or "").strip()
    position_vcf = (row.get("PositionVCF") or "").strip()
    alleles = (row.get("ReferenceAlleleVCF") or "", row.get("AlternateAlleleVCF") or "")
    # ClinVar marks missing VCF coordinates with -1 and missing alleles with "na".
    if (
        chromosome not in CHROMOSOMES
        or not position_vcf.isdigit()
        or int(position_vcf) < 1
        or not all(len(a.strip()) == 1 and a.strip() in "ACGT" for a in alleles)
    ):
        return None, "missing_coordinates"

    ref_aa, position, alt_aa = change
    rs = (row.get("RS# (dbSNP)") or "").strip()
    record = {
        "VariationID": row["VariationID"].strip(),
        "AlleleID": row["#AlleleID"].strip(),
        "gene": gene,
        "GeneID": row.get("GeneID", "").strip(),
        "HGNC_ID": row.get("HGNC_ID", "").strip(),
        "transcript": parse_transcript(name),
        "protein_change": f"{ref_aa}{position}{alt_aa}",
        "ref_aa": ref_aa,
        "protein_position": position,
        "alt_aa": alt_aa,
        "label": LABELS[significance],
        "ClinicalSignificance": significance,
        "ReviewStatus": review,
        "review_stars": stars,
        "NumberSubmitters": int(row.get("NumberSubmitters") or 0),
        "LastEvaluated": row.get("LastEvaluated", "").strip(),
        "Assembly": row["Assembly"],
        "Chromosome": chromosome,
        "Start": row.get("Start", "").strip(),
        "Stop": row.get("Stop", "").strip(),
        "PositionVCF": row["PositionVCF"].strip(),
        "ReferenceAlleleVCF": row["ReferenceAlleleVCF"].strip(),
        "AlternateAlleleVCF": row["AlternateAlleleVCF"].strip(),
        "rsID": f"rs{rs}" if usable(rs) and rs != "-1" else "",
        "Name": name,
    }
    return record, ""


def scan_clinvar(rows: Iterable[dict[str, str]]) -> tuple[list[dict[str, Any]], Counter]:
    """Filter ClinVar rows. Returns kept records and a Counter of drop reasons."""
    kept: list[dict[str, Any]] = []
    reasons: Counter = Counter()
    for row in rows:
        reasons["rows_read"] += 1
        record, reason = filter_row(row)
        if record is None:
            reasons[reason] += 1
        else:
            kept.append(record)
    reasons["rows_kept"] = len(kept)
    return kept, reasons


def read_clinvar(path: Path) -> Iterator[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _quality_key(record: dict[str, Any]) -> tuple:
    """Sort key preferring better-supported records; VariationID breaks ties."""
    return (-record["review_stars"], -record["NumberSubmitters"], int(record["VariationID"]))


def deduplicate(records: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], Counter]:
    """Remove duplicate VariationIDs and duplicate protein changes within a gene.

    Distinct nucleotide changes can give the same amino-acid change. Only the
    best-supported record is kept, and if the records disagree on the label
    the protein change is dropped altogether.
    """
    stats: Counter = Counter()
    by_variation: dict[str, dict[str, Any]] = {}
    for record in sorted(records, key=_quality_key):
        if record["VariationID"] in by_variation:
            stats["duplicate_variation_id"] += 1
            continue
        by_variation[record["VariationID"]] = record

    by_change: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for record in by_variation.values():
        by_change[(record["gene"], record["protein_change"])].append(record)

    unique: list[dict[str, Any]] = []
    for group in by_change.values():
        if len({r["label"] for r in group}) > 1:
            stats["protein_change_label_conflict"] += len(group)
            continue
        group.sort(key=_quality_key)
        stats["duplicate_protein_change"] += len(group) - 1
        unique.append(group[0])
    unique.sort(key=lambda r: (r["gene"], r["protein_position"], r["protein_change"]))
    return unique, stats


def class_counts(records: Iterable[dict[str, Any]]) -> dict[str, Counter]:
    counts: dict[str, Counter] = defaultdict(Counter)
    for record in records:
        counts[record["gene"]][record["label"]] += 1
    return counts


def eligible_genes(records: list[dict[str, Any]], per_class: int) -> list[str]:
    counts = class_counts(records)
    return sorted(g for g, c in counts.items() if c[1] >= per_class and c[0] >= per_class)


# ---------------------------------------------------------------------------
# UniProt / AlphaFold mapping
# ---------------------------------------------------------------------------

@dataclass
class ApiClient:
    """HTTP client with an on-disk JSON cache and retries."""

    cache_dir: Path
    offline: bool = False
    delay: float = 0.1
    session: requests.Session = field(default_factory=requests.Session)

    def get_json(self, cache_name: str, url: str, params: dict[str, str] | None = None) -> Any:
        cache_path = self.cache_dir / cache_name
        if cache_path.exists():
            return json.loads(cache_path.read_text(encoding="utf-8"))
        if self.offline:
            raise RuntimeError(f"offline mode and no cache entry: {cache_path}")
        for attempt in range(4):
            try:
                response = self.session.get(url, params=params, timeout=30)
                if response.status_code == 404:
                    data: Any = None
                    break
                response.raise_for_status()
                data = response.json()
                break
            except (requests.RequestException, ValueError):
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
        time.sleep(self.delay)
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(data) + "\n", encoding="utf-8")
        return data


def uniprot_entry(client: ApiClient, gene: str) -> tuple[dict[str, Any] | None, str]:
    """The single reviewed human UniProt entry whose primary gene name is ``gene``."""
    data = client.get_json(
        f"uniprot_{gene}.json",
        UNIPROT_SEARCH,
        {
            "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
            "fields": "accession,gene_primary,length,sequence",
            "format": "json",
            "size": "25",
        },
    )
    matches = [
        r for r in (data or {}).get("results", [])
        if any(g.get("geneName", {}).get("value") == gene for g in r.get("genes", [])[:1])
    ]
    if not matches:
        return None, "no_reviewed_uniprot_entry"
    if len(matches) > 1:
        return None, "multiple_uniprot_entries"
    entry = matches[0]
    return {
        "accession": entry["primaryAccession"],
        "sequence": entry["sequence"]["value"],
        "length": entry["sequence"]["length"],
    }, ""


def alphafold_model(client: ApiClient, accession: str, sequence: str) -> tuple[dict[str, Any] | None, str]:
    """The full-length canonical AlphaFold DB model for ``accession``."""
    data = client.get_json(f"alphafold_{accession}.json", ALPHAFOLD_API.format(accession=accession))
    if not data:
        return None, "no_alphafold_model"
    models = [m for m in data if m.get("uniprotAccession") == accession]
    if not models:
        return None, "no_canonical_alphafold_model"
    model = models[0]
    model_sequence = model.get("uniprotSequence") or model.get("sequence") or ""
    if model.get("uniprotStart") != 1 or model.get("uniprotEnd") != len(sequence):
        return None, "alphafold_model_not_full_length"
    if model_sequence != sequence:
        return None, "alphafold_sequence_differs_from_uniprot"
    return {
        "alphafold_entry": model["entryId"],
        "alphafold_version": model.get("latestVersion"),
        "alphafold_mean_plddt": model.get("globalMetricValue"),
        "alphafold_cif_url": model.get("cifUrl"),
        "alphafold_sequence_checksum": model.get("sequenceChecksum", ""),
    }, ""


def map_gene(client: ApiClient, gene: str, max_length: int) -> tuple[dict[str, Any] | None, str]:
    entry, reason = uniprot_entry(client, gene)
    if entry is None:
        return None, reason
    if entry["length"] > max_length:
        return None, "protein_longer_than_max_length"
    model, reason = alphafold_model(client, entry["accession"], entry["sequence"])
    if model is None:
        return None, reason
    return {**entry, **model}, ""


def residue_matches(record: dict[str, Any], sequence: str) -> bool:
    position = record["protein_position"]
    return 1 <= position <= len(sequence) and sequence[position - 1] == record["ref_aa"]


# ---------------------------------------------------------------------------
# Gene ranking and sampling
# ---------------------------------------------------------------------------

def rank_genes(counts: dict[str, Counter], genes: Iterable[str]) -> list[dict[str, Any]]:
    """Rank genes by the size of their smaller class, then class balance, then name.

    Genes with more variants in their minority class give the random sample
    more to choose from, which makes the per-gene sample less dependent on
    a handful of records.
    """
    rows = []
    for gene in genes:
        pathogenic, benign = counts[gene][1], counts[gene][0]
        rows.append({
            "gene": gene,
            "available_pathogenic": pathogenic,
            "available_benign": benign,
            "minority_class_count": min(pathogenic, benign),
            "class_balance": round(min(pathogenic, benign) / max(pathogenic, benign), 4),
        })
    rows.sort(key=lambda r: (-r["minority_class_count"], -r["class_balance"], r["gene"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows


def sample_gene(records: list[dict[str, Any]], gene: str, label: int, k: int, seed: int) -> list[dict[str, Any]]:
    """Seeded random sample of ``k`` records for one gene and class.

    Each (gene, label) pair has its own RNG so that changing the gene list
    does not change the variants drawn for any other gene.
    """
    pool = sorted(
        (r for r in records if r["gene"] == gene and r["label"] == label),
        key=lambda r: int(r["VariationID"]),
    )
    rng = random.Random(f"{seed}:{gene}:{label}")
    return sorted(rng.sample(pool, k), key=lambda r: int(r["VariationID"]))


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate_cohort(
    cohort: list[dict[str, Any]],
    genes: list[dict[str, Any]],
    sequences: dict[str, str],
    target_genes: int,
    per_class: int,
) -> list[dict[str, Any]]:
    """Run every acceptance check. Returns a list of {check, passed, detail}."""
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: Any = "") -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    gene_names = [g["gene"] for g in genes]
    counts = class_counts(cohort)
    bad_genes = {g: dict(c) for g, c in counts.items() if c[1] != per_class or c[0] != per_class}
    labels = Counter(r["label"] for r in cohort)

    check("gene_count", len(gene_names) == target_genes, len(gene_names))
    check("genes_distinct", len(set(gene_names)) == len(gene_names))
    check("cohort_genes_match_selected_genes", set(counts) == set(gene_names))
    check("per_gene_class_counts", not bad_genes, bad_genes or f"{per_class} pathogenic + {per_class} benign each")
    check("total_variants", len(cohort) == target_genes * per_class * 2, len(cohort))
    check("class_balance", labels[1] == labels[0], {"pathogenic": labels[1], "benign": labels[0]})

    def unique(key_fn) -> bool:
        keys = [key_fn(r) for r in cohort]
        return len(keys) == len(set(keys))

    check("unique_variant_id", unique(lambda r: r["variant_id"]))
    check("unique_variation_id", unique(lambda r: r["VariationID"]))
    check("unique_genomic_allele", unique(lambda r: (r["Chromosome"], r["PositionVCF"], r["ReferenceAlleleVCF"], r["AlternateAlleleVCF"])))
    check("unique_protein_change_per_gene", unique(lambda r: (r["gene"], r["protein_change"])))
    check("assembly_grch38_only", all(r["Assembly"] == config.ASSEMBLY for r in cohort))
    check("missense_only", all(parse_protein_change(r["Name"]) == (r["ref_aa"], r["protein_position"], r["alt_aa"]) for r in cohort))
    check("labels_consistent", all(LABELS[r["ClinicalSignificance"]] == r["label"] for r in cohort))
    check("review_stars_at_least_one", all(r["review_stars"] >= MIN_REVIEW_STARS for r in cohort))
    check("uniprot_accessions_distinct", len({g["accession"] for g in genes}) == len(genes))
    check("alphafold_model_for_every_gene", all(g.get("alphafold_entry") for g in genes))
    check("protein_length_within_limit", all(g["length"] <= config.MAX_PROTEIN_LENGTH for g in genes))
    accession = {g["gene"]: g["accession"] for g in genes}
    mismatched = [r["variant_id"] for r in cohort if not residue_matches(r, sequences[accession[r["gene"]]])]
    check("reference_residue_matches_alphafold_sequence", not mismatched, mismatched[:10])
    return checks


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def _git_commit() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=config.PROJECT_ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def _write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_cohort(
    clinvar_path: Path = config.CLINVAR_FILE,
    output_dir: Path = config.PHASE1_DIR,
    cache_dir: Path = config.PHASE1_CACHE_DIR,
    target_genes: int = config.TARGET_GENES,
    per_class: int = config.VARIANTS_PER_CLASS,
    seed: int = config.SEED,
    offline: bool = False,
) -> dict[str, Any]:
    if not clinvar_path.exists():
        raise FileNotFoundError(f"ClinVar file not found: {clinvar_path}. Download it from {config.CLINVAR_URL}")

    print(f"Scanning {clinvar_path} ...", flush=True)
    records, scan_stats = scan_clinvar(read_clinvar(clinvar_path))
    print(f"  kept {len(records):,} of {scan_stats['rows_read']:,} rows", flush=True)

    records, dedup_stats = deduplicate(records)
    print(f"  {len(records):,} unique missense variants after deduplication", flush=True)

    candidates = eligible_genes(records, per_class)
    print(f"  {len(candidates)} genes with >= {per_class} variants in both classes", flush=True)

    client = ApiClient(cache_dir=cache_dir, offline=offline)
    mapping: dict[str, dict[str, Any]] = {}
    candidate_rows: list[dict[str, Any]] = []
    residue_stats: Counter = Counter()
    kept_records: list[dict[str, Any]] = []
    by_gene: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for record in records:
        by_gene[record["gene"]].append(record)

    for index, gene in enumerate(candidates, start=1):
        if index % 50 == 0:
            print(f"  mapped {index}/{len(candidates)} genes", flush=True)
        info, reason = map_gene(client, gene, config.MAX_PROTEIN_LENGTH)
        row = {"gene": gene, "mapping_status": "mapped" if info else "excluded", "mapping_reason": reason}
        before = Counter(r["label"] for r in by_gene[gene])
        row.update({"clinvar_pathogenic": before[1], "clinvar_benign": before[0]})
        if info:
            matched = [r for r in by_gene[gene] if residue_matches(r, info["sequence"])]
            residue_stats["residue_match"] += len(matched)
            residue_stats["residue_mismatch"] += len(by_gene[gene]) - len(matched)
            after = Counter(r["label"] for r in matched)
            row.update({
                "uniprot_accession": info["accession"],
                "protein_length": info["length"],
                "residue_matched_pathogenic": after[1],
                "residue_matched_benign": after[0],
            })
            if after[1] < per_class or after[0] < per_class:
                row.update(mapping_status="excluded", mapping_reason="too_few_variants_after_residue_check")
            else:
                mapping[gene] = info
                kept_records.extend(matched)
        candidate_rows.append(row)

    counts = class_counts(kept_records)
    ranked = rank_genes(counts, mapping)
    print(f"  {len(ranked)} genes pass structure mapping", flush=True)
    if len(ranked) < target_genes:
        raise RuntimeError(f"only {len(ranked)} eligible genes, need {target_genes}")
    selected = ranked[:target_genes]

    for row in candidate_rows:
        row["selection_rank"] = next((r["rank"] for r in ranked if r["gene"] == row["gene"]), "")
        row["selected"] = row["gene"] in {s["gene"] for s in selected}

    cohort: list[dict[str, Any]] = []
    gene_rows: list[dict[str, Any]] = []
    for item in selected:
        gene = item["gene"]
        info = mapping[gene]
        for label in (1, 0):
            for record in sample_gene(kept_records, gene, label, per_class, seed):
                cohort.append({
                    **record,
                    "variant_id": f"{gene}:p.{record['protein_change']}",
                    "uniprot_accession": info["accession"],
                })
        gene_rows.append({
            **item,
            "uniprot_accession": info["accession"],
            "protein_length": info["length"],
            "alphafold_entry": info["alphafold_entry"],
            "alphafold_version": info["alphafold_version"],
            "alphafold_mean_plddt": info["alphafold_mean_plddt"],
            "alphafold_cif_url": info["alphafold_cif_url"],
            "selected_pathogenic": per_class,
            "selected_benign": per_class,
        })
    cohort.sort(key=lambda r: (r["gene"], r["protein_position"], r["protein_change"]))

    genes_for_checks = [{**g, "accession": g["uniprot_accession"], "length": g["protein_length"]} for g in gene_rows]
    sequences = {mapping[g["gene"]]["accession"]: mapping[g["gene"]]["sequence"] for g in gene_rows}
    checks = validate_cohort(cohort, genes_for_checks, sequences, target_genes, per_class)
    passed = all(c["passed"] for c in checks)

    output_dir.mkdir(parents=True, exist_ok=True)
    cohort_path = output_dir / "cohort.csv"
    genes_path = output_dir / "selected_genes.csv"
    candidates_path = output_dir / "gene_candidates.csv"
    _write_csv(cohort_path, cohort, COHORT_FIELDS)
    _write_csv(genes_path, gene_rows, [
        "rank", "gene", "uniprot_accession", "protein_length", "alphafold_entry", "alphafold_version",
        "alphafold_mean_plddt", "alphafold_cif_url", "available_pathogenic", "available_benign",
        "minority_class_count", "class_balance", "selected_pathogenic", "selected_benign",
    ])
    _write_csv(candidates_path, candidate_rows, [
        "gene", "mapping_status", "mapping_reason", "uniprot_accession", "protein_length",
        "clinvar_pathogenic", "clinvar_benign", "residue_matched_pathogenic", "residue_matched_benign",
        "selection_rank", "selected",
    ])

    lengths = sorted(g["protein_length"] for g in gene_rows)
    metadata = {
        "phase": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": _git_commit(),
        "source": {
            "clinvar_url": config.CLINVAR_URL,
            "clinvar_file": clinvar_path.relative_to(config.PROJECT_ROOT).as_posix()
            if clinvar_path.is_relative_to(config.PROJECT_ROOT) else str(clinvar_path),
            "clinvar_sha256": sha256_file(clinvar_path),
            "clinvar_file_modified_utc": datetime.fromtimestamp(clinvar_path.stat().st_mtime, timezone.utc).isoformat(timespec="seconds"),
            "uniprot_api": UNIPROT_SEARCH,
            "alphafold_api": ALPHAFOLD_API,
        },
        "parameters": {
            "target_genes": target_genes,
            "variants_per_class": per_class,
            "seed": seed,
            "assembly": config.ASSEMBLY,
            "variant_type": "single nucleotide variant",
            "origin": "germline",
            "labels": LABELS,
            "review_status_stars": REVIEW_STARS,
            "min_review_stars": MIN_REVIEW_STARS,
            "max_protein_length": config.MAX_PROTEIN_LENGTH,
            "gene_ranking": "minority-class count desc, class balance desc, gene asc",
            "sampling": "random.Random(f'{seed}:{gene}:{label}').sample over VariationID-sorted pool",
        },
        "filter_counts": dict(scan_stats),
        "deduplication": dict(dedup_stats),
        "unique_missense_variants": len(records),
        "genes_with_enough_variants": len(candidates),
        "gene_mapping_exclusions": dict(Counter(r["mapping_reason"] for r in candidate_rows if r["mapping_status"] == "excluded")),
        "residue_check": dict(residue_stats),
        "genes_passing_mapping": len(ranked),
        "cohort": {
            "variants": len(cohort),
            "genes": len(gene_rows),
            "pathogenic": sum(r["label"] == 1 for r in cohort),
            "benign": sum(r["label"] == 0 for r in cohort),
            "review_stars": dict(sorted(Counter(r["review_stars"] for r in cohort).items())),
            "clinical_significance": dict(Counter(r["ClinicalSignificance"] for r in cohort)),
            "protein_length_min": lengths[0],
            "protein_length_median": lengths[len(lengths) // 2],
            "protein_length_max": lengths[-1],
            "cohort_sha256": sha256_file(cohort_path),
        },
    }
    write_json(output_dir / "metadata.json", metadata)
    write_json(output_dir / "validation_report.json", {"passed": passed, "checks": checks})

    print(f"\nCohort: {len(cohort):,} variants, {len(gene_rows)} genes")
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['check']}")
    return {"passed": passed, "metadata": metadata, "checks": checks}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--clinvar", type=Path, default=config.CLINVAR_FILE)
    parser.add_argument("--output-dir", type=Path, default=config.PHASE1_DIR)
    parser.add_argument("--cache-dir", type=Path, default=config.PHASE1_CACHE_DIR)
    parser.add_argument("--target-genes", type=int, default=config.TARGET_GENES)
    parser.add_argument("--variants-per-class", type=int, default=config.VARIANTS_PER_CLASS)
    parser.add_argument("--seed", type=int, default=config.SEED)
    parser.add_argument("--offline", action="store_true", help="use only cached API responses")
    args = parser.parse_args(argv)
    result = build_cohort(
        clinvar_path=args.clinvar,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        target_genes=args.target_genes,
        per_class=args.variants_per_class,
        seed=args.seed,
        offline=args.offline,
    )
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
