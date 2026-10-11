import numpy as np
import pandas as pd

from biodb.phase7_stability import gene_bootstrap_auc_difference


def _oof(probabilities):
    rng = np.random.default_rng(1)
    labels = np.tile([0, 1], 100)
    return pd.DataFrame({"gene": np.repeat([f"G{i}" for i in range(20)], 10), "label": labels,
                         "m": probabilities(labels, rng)})


def test_bootstrap_detects_clear_improvement():
    weak = _oof(lambda y, rng: y * 0.2 + rng.normal(0, 0.5, len(y)))
    strong = weak.assign(m=weak["label"] * 2 + np.random.default_rng(2).normal(0, 0.5, len(weak)))
    result = gene_bootstrap_auc_difference(weak, strong, "m")
    assert result["delta_roc_auc"] > 0.2
    assert result["ci95_low"] > 0


def test_bootstrap_of_identical_predictions_is_zero():
    oof = _oof(lambda y, rng: y + rng.normal(0, 1, len(y)))
    result = gene_bootstrap_auc_difference(oof, oof.copy(), "m")
    assert result["delta_roc_auc"] == 0 and result["ci95_low"] == 0 and result["ci95_high"] == 0
