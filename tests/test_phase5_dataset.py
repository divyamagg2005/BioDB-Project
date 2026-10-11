import csv
import json

import pytest

from biodb import config
from biodb.phase5_dataset import ALL_FEATURES, FEATURE_DESCRIPTIONS, build_row, substitution_features

COHORT = {"variant_id": "TP53:p.Y220C", "gene": "TP53", "uniprot_accession": "P04637", "protein_position": "220",
          "ref_aa": "Y", "alt_aa": "C", "label": "1"}


def test_substitution_features_y220c():
    f = substitution_features("Y", "C")
    assert f["blosum62"] == -2
    assert f["hydropathy_change"] == pytest.approx(3.8)
    assert f["volume_change"] == pytest.approx(-85.1)
    assert f["to_cysteine"] == 1 and f["from_cysteine"] == 0
    assert substitution_features("D", "K")["charge_change"] == 2


def test_build_row_keeps_missing_and_flags_gnomad_absence():
    row = build_row(COHORT, {"revel": "", "cadd_phred": "25.1"}, {"plddt": "90", "secondary_structure": "strand"},
                    {"ddg": "3.48", "status": "ok"})
    assert row["revel"] is None
    assert row["cadd_phred"] == 25.1
    assert row["gnomad_absent"] == 1 and row["gnomad_af_max"] is None
    assert row["ss_strand"] == 1 and row["ss_helix"] == 0
    assert row["ddg"] == 3.48
    assert row["label"] == 1


def test_failed_foldx_gives_missing_ddg():
    row = build_row(COHORT, {}, {}, {"ddg": "1.0", "status": "error"})
    assert row["ddg"] is None


def test_every_feature_documented():
    assert set(ALL_FEATURES) == set(FEATURE_DESCRIPTIONS)


VALIDATION = config.PHASE5_DIR / "validation_report.json"


@pytest.mark.data
def test_output_validation_passed():
    if not VALIDATION.exists():
        pytest.skip("Phase 5 outputs not generated yet")
    report = json.loads(VALIDATION.read_text())
    assert report["passed"], [c for c in report["checks"] if not c["passed"]]
