import pandas as pd

from biodb.phase8_analysis import auc


def test_auc_ignores_missing_and_needs_both_classes():
    y = pd.Series([0, 0, 1, 1, 1])
    assert auc(y, pd.Series([0.1, 0.2, 0.8, 0.9, None])) == 1.0
    assert auc(pd.Series([1, 1]), pd.Series([0.3, 0.4])) is None
