import csv
import json
from pathlib import Path

import pytest

from biodb import config
from biodb.phase3_structure import residue_features, secondary_structure

# Minimal single-chain mmCIF: three residues (ALA, GLY, ARG) with pLDDT in B_iso
# and one helix record covering residues 1-2.
MINI_CIF = """data_test
#
loop_
_struct_conf.conf_type_id
_struct_conf.beg_label_seq_id
_struct_conf.end_label_seq_id
HELX_RH_AL_P 1 2
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 1 ? 0.000 0.000 0.000 1.00 90.00 ? 1 ALA A N 1
ATOM 2 C CA . ALA A 1 1 ? 1.458 0.000 0.000 1.00 90.00 ? 1 ALA A CA 1
ATOM 3 C C . ALA A 1 1 ? 2.009 1.420 0.000 1.00 90.00 ? 1 ALA A C 1
ATOM 4 O O . ALA A 1 1 ? 1.251 2.390 0.000 1.00 90.00 ? 1 ALA A O 1
ATOM 5 C CB . ALA A 1 1 ? 1.988 -0.773 -1.199 1.00 90.00 ? 1 ALA A CB 1
ATOM 6 N N . GLY A 1 2 ? 3.332 1.536 0.000 1.00 60.00 ? 2 GLY A N 1
ATOM 7 C CA . GLY A 1 2 ? 3.988 2.839 0.000 1.00 60.00 ? 2 GLY A CA 1
ATOM 8 C C . GLY A 1 2 ? 5.504 2.693 0.000 1.00 60.00 ? 2 GLY A C 1
ATOM 9 O O . GLY A 1 2 ? 6.030 1.580 0.000 1.00 60.00 ? 2 GLY A O 1
ATOM 10 N N . ARG A 1 3 ? 6.210 3.813 0.000 1.00 30.00 ? 3 ARG A N 1
ATOM 11 C CA . ARG A 1 3 ? 7.665 3.813 0.000 1.00 30.00 ? 3 ARG A CA 1
ATOM 12 C C . ARG A 1 3 ? 8.216 5.233 0.000 1.00 30.00 ? 3 ARG A C 1
ATOM 13 O O . ARG A 1 3 ? 7.458 6.203 0.000 1.00 30.00 ? 3 ARG A O 1
ATOM 14 C CB . ARG A 1 3 ? 8.195 3.040 -1.199 1.00 30.00 ? 3 ARG A CB 1
#
"""


@pytest.fixture
def mini_cif(tmp_path: Path) -> Path:
    path = tmp_path / "mini.cif"
    path.write_text(MINI_CIF)
    return path


def test_secondary_structure_maps_dssp_records():
    cif = {
        "_struct_conf.conf_type_id": ["HELX_RH_AL_P", "STRN", "BEND"],
        "_struct_conf.beg_label_seq_id": ["1", "4", "6"],
        "_struct_conf.end_label_seq_id": ["2", "5", "6"],
    }
    assert secondary_structure(cif, 6) == ["helix", "helix", "coil", "strand", "strand", "coil"]


def test_residue_features_on_minimal_model(mini_cif):
    sequence, features = residue_features(mini_cif)
    assert sequence == "AGR"
    assert features[1]["plddt"] == 90.0
    assert features[2]["plddt_window"] == 60.0
    assert features[1]["secondary_structure"] == "helix"
    assert features[3]["secondary_structure"] == "coil"
    # Glycine has no CB, so the CA is used; all three residues are within 10 A.
    assert all(f["neighbor_count"] == 2 for f in features.values())
    assert all(0 < f["rsa"] <= 1.5 for f in features.values())


VALIDATION = config.PHASE3_DIR / "validation_report.json"
FEATURES = config.PHASE3_DIR / "variant_structure_features.csv"


@pytest.mark.data
def test_output_validation_passed():
    if not VALIDATION.exists():
        pytest.skip("Phase 3 outputs not generated yet")
    report = json.loads(VALIDATION.read_text())
    assert report["passed"], [c for c in report["checks"] if not c["passed"]]


@pytest.mark.data
def test_output_wt_residues_match_cohort():
    if not FEATURES.exists():
        pytest.skip("Phase 3 outputs not generated yet")
    with FEATURES.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with (config.PHASE1_DIR / "cohort.csv").open(encoding="utf-8") as handle:
        cohort = {r["variant_id"]: r for r in csv.DictReader(handle)}
    assert len(rows) == len(cohort)
    assert all(r["wt_aa"] == cohort[r["variant_id"]]["ref_aa"] for r in rows)
