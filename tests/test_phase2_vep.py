import csv
import json

import pytest

from biodb import config
from biodb.phase2_vep import allele_frequencies, choose_transcript, parse_result, vep_input

ROW = {
    "variant_id": "TP53:p.R175H", "gene": "TP53", "Chromosome": "17", "PositionVCF": "7675088",
    "ReferenceAlleleVCF": "C", "AlternateAlleleVCF": "T", "ref_aa": "R", "alt_aa": "H",
    "protein_position": "175",
}


def transcript(**overrides):
    t = {
        "gene_symbol": "TP53", "transcript_id": "ENST1", "protein_start": 175, "amino_acids": "R/H",
        "consequence_terms": ["missense_variant"], "cadd_phred": 29.1, "revel": 0.93,
        "alphamissense": {"am_pathogenicity": 0.99, "am_class": "likely_pathogenic"},
        "sift_score": 0.0, "polyphen_score": 1.0, "conservation": 4.5,
    }
    t.update(overrides)
    return t


def result(consequences, **overrides):
    r = {"id": ROW["variant_id"], "seq_region_name": "17", "start": 7675088, "allele_string": "C/T",
         "transcript_consequences": consequences}
    r.update(overrides)
    return r


def test_vep_input_uses_vcf_position_and_variant_id():
    assert vep_input(ROW) == "17 7675088 TP53:p.R175H C T . . ."


def test_choose_transcript_prefers_mane_then_canonical():
    other = transcript(transcript_id="ENST3")
    canonical = transcript(transcript_id="ENST2", canonical=1)
    mane = transcript(transcript_id="ENST9", mane_select="NM_000546.6")
    assert choose_transcript([other, canonical, mane], ROW) == (mane, "mane_select")
    assert choose_transcript([other, canonical], ROW) == (canonical, "canonical")
    assert choose_transcript([other], ROW) == (other, "other_matching")


def test_choose_transcript_requires_matching_protein_change():
    wrong_position = transcript(protein_start=174, canonical=1)
    wrong_gene = transcript(gene_symbol="WRAP53")
    assert choose_transcript([wrong_position], ROW) == (None, "no_transcript_matches_protein_change")
    assert choose_transcript([wrong_gene], ROW) == (None, "no_transcript_for_gene")


def test_allele_frequencies_use_exact_alt_allele():
    r = {"colocated_variants": [{"frequencies": {"T": {"gnomade": 1e-5, "gnomadg": 3e-5}, "G": {"gnomade": 0.2}}}]}
    assert allele_frequencies(r, "T") == {"gnomade_af": 1e-5, "gnomadg_af": 3e-5, "gnomad_af_max": 3e-5}
    assert allele_frequencies({}, "T") == {"gnomade_af": None, "gnomadg_af": None, "gnomad_af_max": None}


def test_parse_result_extracts_features_and_ignores_clinvar():
    r = result([transcript(mane_select="NM_000546.6")],
               colocated_variants=[{"clin_sig": ["pathogenic"], "frequencies": {"T": {"gnomadg": 1e-6}}}])
    parsed = parse_result(r, ROW)
    assert parsed["annotation_status"] == "annotated"
    assert parsed["revel"] == 0.93 and parsed["alphamissense"] == 0.99
    assert parsed["gnomadg_af"] == 1e-6 and parsed["gnomade_af"] is None
    assert "pathogenic" not in json.dumps(parsed).replace("likely_pathogenic", "")


def test_parse_result_rejects_coordinate_and_allele_mismatch():
    assert parse_result(result([transcript()], start=1), ROW)["failure_reason"] == "coordinate_mismatch"
    assert parse_result(result([transcript()], allele_string="C/A"), ROW)["failure_reason"] == "allele_mismatch"
    assert parse_result(None, ROW)["failure_reason"] == "no_vep_result"


def test_parse_result_keeps_missing_scores_missing():
    parsed = parse_result(result([transcript(revel=None, alphamissense=None)]), ROW)
    assert parsed["revel"] is None and parsed["alphamissense"] is None


VALIDATION = config.PHASE2_DIR / "validation_report.json"
ANNOTATIONS = config.PHASE2_DIR / "vep_annotations.csv"


@pytest.mark.data
def test_output_validation_passed():
    if not VALIDATION.exists():
        pytest.skip("Phase 2 outputs not generated yet")
    report = json.loads(VALIDATION.read_text())
    assert report["passed"], [c for c in report["checks"] if not c["passed"]]


@pytest.mark.data
def test_output_has_one_row_per_cohort_variant():
    if not ANNOTATIONS.exists():
        pytest.skip("Phase 2 outputs not generated yet")
    with ANNOTATIONS.open(encoding="utf-8") as handle:
        ids = [r["variant_id"] for r in csv.DictReader(handle)]
    with (config.PHASE1_DIR / "cohort.csv").open(encoding="utf-8") as handle:
        cohort_ids = [r["variant_id"] for r in csv.DictReader(handle)]
    assert sorted(ids) == sorted(cohort_ids)


def test_batch_lines_sorted_by_chromosome_then_position():
    rows = [
        {**ROW, "variant_id": "a", "Chromosome": "X", "PositionVCF": "5"},
        {**ROW, "variant_id": "b", "Chromosome": "17", "PositionVCF": "20"},
        {**ROW, "variant_id": "c", "Chromosome": "2", "PositionVCF": "9"},
        {**ROW, "variant_id": "d", "Chromosome": "17", "PositionVCF": "3"},
    ]
    from biodb.phase2_vep import batch_lines
    assert [line.split()[2] for line in batch_lines(rows)] == ["c", "d", "b", "a"]


def test_missing_score_filled_from_transcript_with_same_change():
    chosen = transcript(mane_select="NM_000546.6", revel=None, polyphen_score=None)
    peer_a = transcript(transcript_id="ENST2", revel=0.8, polyphen_score=0.9)
    # Isoforms may number the residue differently, so position is not compared.
    peer_b = transcript(transcript_id="ENST3", revel=0.6, polyphen_score=None, protein_start=42)
    unrelated = transcript(transcript_id="ENST4", amino_acids="R/C", revel=0.1, polyphen_score=0.1)
    parsed = parse_result(result([chosen, peer_a, peer_b, unrelated]), ROW)
    assert parsed["revel"] == pytest.approx(0.7)
    assert parsed["polyphen_score"] == 0.9
    assert parsed["scores_from_other_transcript"] == "revel;polyphen_score"
    assert parse_result(result([transcript(mane_select="x")]), ROW)["scores_from_other_transcript"] == ""
