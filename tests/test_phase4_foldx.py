import csv
import json
import shutil
from pathlib import Path

import pytest

from biodb import config
from biodb.phase4_foldx import foldx_mutation, parse_buildmodel, read_fxout_table

# Real FoldX 5.1 BuildModel output for TP53 Y220C (3 runs) on the repaired AlphaFold model.
FIXTURES = Path(__file__).parent / "fixtures" / "foldx"


def test_foldx_mutation_includes_chain():
    assert foldx_mutation("Y", 220, "C") == "YA220C"


def test_read_fxout_table_skips_banner():
    rows = read_fxout_table(FIXTURES / "Average_TP53_Repair.fxout")
    assert len(rows) == 1
    assert rows[0]["Pdb"] == "TP53_Repair_1"
    assert rows[0]["total energy"] == "3.47969"


def test_parse_buildmodel_reads_total_energy_not_last_column():
    result = parse_buildmodel(FIXTURES, "TP53_Repair")
    # Total-energy ΔΔG, not the trailing "Entropy Complex" column (0).
    assert result["ddg"] == pytest.approx(3.47969)
    assert result["ddg_sd"] == pytest.approx(0.0371959)
    assert result["wt_energy"] == pytest.approx(214.3843, abs=1e-3)
    assert result["mutant_energy"] == pytest.approx(217.8643, abs=1e-3)
    assert result["mutant_energy"] - result["wt_energy"] == pytest.approx(result["ddg"], abs=0.01)


def test_parse_buildmodel_rejects_missing_output(tmp_path):
    shutil.copy(FIXTURES / "Average_TP53_Repair.fxout", tmp_path)
    with pytest.raises(FileNotFoundError):
        parse_buildmodel(tmp_path, "TP53_Repair")


DDG = config.PHASE4_DIR / "foldx_ddg.csv"


@pytest.mark.data
def test_output_has_one_row_per_cohort_variant():
    if not DDG.exists():
        pytest.skip("Phase 4 outputs not generated yet")
    with DDG.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with (config.PHASE1_DIR / "cohort.csv").open(encoding="utf-8") as handle:
        cohort = [r["variant_id"] for r in csv.DictReader(handle)]
    assert sorted(r["variant_id"] for r in rows) == sorted(cohort)
    assert all(r["status"] in {"ok", "error", "missing"} for r in rows)
    assert all(r["error"] for r in rows if r["status"] != "ok")
