import numpy as np
import pandas as pd
import pytest

from biodb import modeling


@pytest.fixture
def toy(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    rows = []
    for g in range(20):
        for i in range(10):
            label = i % 2
            rows.append({"variant_id": f"G{g}:{i}", "gene": f"G{g}", "label": label,
                         "signal": label + rng.normal(0, 0.5), "noise": rng.normal(), "gappy": np.nan if i == 3 else 1.0})
    monkeypatch.setattr(modeling, "FOLDS", tmp_path / "folds.csv")
    return pd.DataFrame(rows)


def test_folds_are_gene_grouped_balanced_and_reused(toy):
    folds = modeling.get_folds(toy)
    assert (toy.assign(fold=folds).groupby("gene")["fold"].nunique() == 1).all()
    assert folds.value_counts().tolist() == [40] * 5
    assert (modeling.get_folds(toy) == folds).all()


def test_saved_folds_must_match_dataset(toy):
    modeling.get_folds(toy)
    with pytest.raises(ValueError):
        modeling.get_folds(toy.iloc[:-10])


@pytest.mark.parametrize("model_name", list(modeling.MODEL_PARAMS))
def test_cross_validate_learns_signal_and_handles_missing(toy, model_name):
    result = modeling.cross_validate(toy, ["signal", "noise", "gappy"], model_name, modeling.get_folds(toy))
    assert result["summary"]["roc_auc"]["mean"] > 0.8
    assert len(result["per_fold"]) == modeling.N_FOLDS
    assert sum(result["confusion_matrix"].values()) == len(toy)


def test_metrics_perfect_prediction():
    m = modeling.metrics(np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9]))
    assert m["roc_auc"] == 1.0 and m["f1"] == 1.0 and m["mcc"] == 1.0
