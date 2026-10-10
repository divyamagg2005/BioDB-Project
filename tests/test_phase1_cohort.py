import csv
import json
from collections import Counter
from pathlib import Path

import pytest

from biodb import config
from biodb.phase1_cohort import (
    ApiClient,
    alphafold_model,
    deduplicate,
    eligible_genes,
    filter_row,
    rank_genes,
    residue_matches,
    sample_gene,
    class_counts,
    uniprot_entry,
)


def clinvar_row(**overrides):
    row = {
        "#AlleleID": "1",
        "Type": "single nucleotide variant",
        "Name": "NM_000546.6(TP53):c.524G>A (p.Arg175His)",
        "GeneID": "7157",
        "GeneSymbol": "TP53",
        "HGNC_ID": "HGNC:11998",
        "ClinicalSignificance": "Pathogenic",
        "OriginSimple": "germline",
        "Assembly": "GRCh38",
        "Chromosome": "17",
        "Start": "7675088",
        "Stop": "7675088",
        "ReviewStatus": "criteria provided, multiple submitters, no conflicts",
        "NumberSubmitters": "5",
        "LastEvaluated": "Jan 01, 2024",
        "RS# (dbSNP)": "28934578",
        "VariationID": "12366",
        "PositionVCF": "7675088",
        "ReferenceAlleleVCF": "C",
        "AlternateAlleleVCF": "T",
    }
    row.update(overrides)
    return row


def test_filter_row_keeps_valid_missense():
    record, reason = filter_row(clinvar_row())
    assert reason == ""
    assert record["protein_change"] == "R175H"
    assert record["label"] == 1
    assert record["review_stars"] == 2
    assert record["rsID"] == "rs28934578"


@pytest.mark.parametrize(
    "overrides, reason",
    [
        ({"Assembly": "GRCh37"}, "other_assembly"),
        ({"ClinicalSignificance": "Uncertain significance"}, "label_not_pathogenic_or_benign"),
        ({"ClinicalSignificance": "Conflicting classifications of pathogenicity"}, "label_not_pathogenic_or_benign"),
        ({"Type": "Deletion"}, "not_snv"),
        ({"OriginSimple": "somatic"}, "not_germline"),
        ({"ReviewStatus": "no assertion criteria provided"}, "review_below_one_star"),
        ({"Name": "NM_000546.6(TP53):c.375G>A (p.Thr125=)"}, "not_missense"),
        ({"Name": "NM_000546.6(TP53):c.586C>T (p.Arg196Ter)"}, "not_missense"),
        ({"Name": "NM_000546.6(TP53):c.2T>A (p.Met1Lys)"}, "start_loss"),
        ({"GeneSymbol": "TP53;WRAP53"}, "ambiguous_gene"),
        ({"Name": "NM_001.1(WRAP53):c.4A>G (p.Ser2Gly)"}, "ambiguous_gene"),
        ({"PositionVCF": "-1"}, "missing_coordinates"),
    ],
)
def test_filter_row_drops(overrides, reason):
    record, got = filter_row(clinvar_row(**overrides))
    assert record is None
    assert got == reason


def record(variation_id, change="R175H", label=1, stars=1, submitters=1, gene="TP53"):
    return {
        "VariationID": str(variation_id),
        "gene": gene,
        "protein_change": change,
        "protein_position": int(change[1:-1]),
        "ref_aa": change[0],
        "label": label,
        "review_stars": stars,
        "NumberSubmitters": submitters,
    }


def test_deduplicate_drops_repeated_variation_id():
    unique, stats = deduplicate([record(1), record(1)])
    assert len(unique) == 1
    assert stats["duplicate_variation_id"] == 1


def test_deduplicate_keeps_best_supported_protein_change():
    unique, stats = deduplicate([record(1, stars=1), record(2, stars=3)])
    assert [r["VariationID"] for r in unique] == ["2"]
    assert stats["duplicate_protein_change"] == 1


def test_deduplicate_drops_label_conflicts():
    unique, stats = deduplicate([record(1, label=1), record(2, label=0), record(3, change="G245S")])
    assert [r["protein_change"] for r in unique] == ["G245S"]
    assert stats["protein_change_label_conflict"] == 2


def test_eligible_genes_needs_both_classes():
    records = [record(i, change=f"A{i}V", label=i % 2) for i in range(1, 7)]
    records += [record(100 + i, change=f"A{i}V", label=1, gene="BRCA1") for i in range(1, 7)]
    assert eligible_genes(records, per_class=3) == ["TP53"]


def test_sample_gene_is_deterministic_and_seed_dependent():
    records = [record(i, change=f"A{i}V", label=1) for i in range(1, 41)]
    first = [r["VariationID"] for r in sample_gene(records, "TP53", 1, 10, seed=42)]
    again = [r["VariationID"] for r in sample_gene(list(reversed(records)), "TP53", 1, 10, seed=42)]
    other = [r["VariationID"] for r in sample_gene(records, "TP53", 1, 10, seed=7)]
    assert first == again
    assert first != other
    assert len(set(first)) == 10


def test_sample_gene_is_not_lowest_coordinates():
    records = [record(i, change=f"A{i}V", label=1) for i in range(1, 101)]
    picked = sorted(int(r["VariationID"]) for r in sample_gene(records, "TP53", 1, 25, seed=42))
    assert picked != list(range(1, 26))


def test_rank_genes_orders_by_minority_class():
    counts = {"A": Counter({1: 30, 0: 30}), "B": Counter({1: 100, 0: 40}), "C": Counter({1: 40, 0: 40})}
    assert [r["gene"] for r in rank_genes(counts, counts)] == ["C", "B", "A"]


def test_residue_matches():
    assert residue_matches(record(1, change="R3H"), "MAR")
    assert not residue_matches(record(1, change="R2H"), "MAR")
    assert not residue_matches(record(1, change="R9H"), "MAR")


def _client_with_cache(tmp_path, name, data):
    (tmp_path / name).write_text(json.dumps(data))
    return ApiClient(cache_dir=tmp_path, offline=True)


def test_uniprot_entry_requires_primary_gene_name(tmp_path):
    data = {"results": [
        {"primaryAccession": "P1", "genes": [{"geneName": {"value": "OTHER"}}], "sequence": {"value": "MA", "length": 2}},
        {"primaryAccession": "P2", "genes": [{"geneName": {"value": "TP53"}}], "sequence": {"value": "MAR", "length": 3}},
    ]}
    entry, reason = uniprot_entry(_client_with_cache(tmp_path, "uniprot_TP53.json", data), "TP53")
    assert reason == ""
    assert entry["accession"] == "P2"


def test_alphafold_model_selects_canonical_full_length(tmp_path):
    canonical = {"uniprotAccession": "P2", "entryId": "AF-P2-F1", "uniprotStart": 1, "uniprotEnd": 3, "uniprotSequence": "MAR"}
    isoform = {**canonical, "uniprotAccession": "P2-2", "entryId": "AF-P2-2-F1"}
    client = _client_with_cache(tmp_path, "alphafold_P2.json", [isoform, canonical])
    model, reason = alphafold_model(client, "P2", "MAR")
    assert reason == "" and model["alphafold_entry"] == "AF-P2-F1"
    model, reason = alphafold_model(client, "P2", "MARK")
    assert model is None and reason == "alphafold_model_not_full_length"


# ---------------------------------------------------------------------------
# Checks on the generated Phase 1 outputs
# ---------------------------------------------------------------------------

COHORT = config.PHASE1_DIR / "cohort.csv"
GENES = config.PHASE1_DIR / "selected_genes.csv"
VALIDATION = config.PHASE1_DIR / "validation_report.json"


def _read(path: Path):
    if not path.exists():
        pytest.skip(f"{path.name} not generated yet")
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


@pytest.mark.data
def test_output_validation_passed():
    if not VALIDATION.exists():
        pytest.skip("validation report not generated yet")
    report = json.loads(VALIDATION.read_text())
    assert report["passed"], [c for c in report["checks"] if not c["passed"]]


@pytest.mark.data
def test_output_cohort_shape():
    cohort = _read(COHORT)
    genes = _read(GENES)
    assert len(genes) == config.TARGET_GENES
    assert len(cohort) == config.TARGET_GENES * config.VARIANTS_PER_CLASS * 2
    counts = class_counts({"gene": r["gene"], "label": int(r["label"])} for r in cohort)
    assert all(c[1] == c[0] == config.VARIANTS_PER_CLASS for c in counts.values())


@pytest.mark.data
def test_output_cohort_has_no_duplicates():
    cohort = _read(COHORT)
    assert len({r["VariationID"] for r in cohort}) == len(cohort)
    assert len({(r["Chromosome"], r["PositionVCF"], r["AlternateAlleleVCF"]) for r in cohort}) == len(cohort)
    assert {r["Assembly"] for r in cohort} == {"GRCh38"}
