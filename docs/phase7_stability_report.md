# Phase 7: Stability-Enhanced Models

## Setup

Run with `python -m biodb.phase7_stability`. Everything is identical to Phase 6: the dataset, the saved gene-grouped folds, the models and hyperparameters, the preprocessing and the seed. Only the feature set changes.

FoldX ΔΔG (`ddg`) is added to four reference feature sets:

| Reference | Features | Question |
|---|---|---|
| `baseline` | conventional + substitution (Phase 6 model) | Does ΔΔG improve the standard model? |
| `baseline_structure` | baseline + AlphaFold structure | Does ΔΔG add anything beyond burial and confidence? |
| `intrinsic` | substitution + structure (no VEP predictors, no gnomAD) | Does ΔΔG help without the circular features? |
| `substitution` | amino-acid physicochemistry only | Does ΔΔG carry signal at all? |

How the comparisons are judged:
- **Primary model:** random forest, chosen in Phase 6 before ΔΔG was added. All three models are reported.
- **Significance:** each comparison gets per-fold paired differences for every metric, plus a 95% CI on the pooled out-of-fold ROC-AUC difference. The CI comes from a cluster bootstrap over genes (1,000 resamples).

## Results: ROC-AUC change from adding ΔΔG

| Reference | Model | ROC-AUC without → with ΔΔG | Δ ROC-AUC [95% CI] | Folds improved | ΔΔG importance rank |
|---|---|---|---|---|---|
| baseline | random forest | 0.978 → 0.978 | +0.0008 [+0.0000, +0.0015] | 4/5 | 8 of 20 |
| baseline | XGBoost | 0.974 → 0.976 | +0.0024 [+0.0011, +0.0037] | 5/5 | 8 of 20 |
| baseline | logistic regression | 0.974 → 0.975 | +0.0008 [−0.0001, +0.0016] | 5/5 | 6 of 20 |
| baseline_structure | random forest | 0.978 → 0.978 | +0.0004 [−0.0001, +0.0009] | 3/5 | 11 of 26 |
| baseline_structure | XGBoost | 0.977 → 0.977 | +0.0001 [−0.0006, +0.0009] | 3/5 | 10 of 26 |
| intrinsic | random forest | 0.886 → 0.888 | +0.0022 [−0.0001, +0.0045] | 5/5 | 2 of 18 |
| intrinsic | XGBoost | 0.886 → 0.886 | +0.0003 [−0.0021, +0.0028] | 4/5 | 5 of 18 |
| substitution | random forest | 0.742 → 0.809 | **+0.0680 [+0.0532, +0.0829]** | 5/5 | 1 of 12 |
| substitution | XGBoost | 0.743 → 0.816 | **+0.0745 [+0.0603, +0.0877]** | 5/5 | 2 of 12 |

Threshold-based metrics follow the same pattern.
- **On the baseline (random forest):** F1 +0.001, MCC +0.002, Brier −0.0013 (better in 5 of 5 folds). Accuracy and recall changes are within fold-to-fold noise.
- **On substitution only (random forest):** accuracy +0.055, precision +0.059, recall +0.056, F1 +0.057, MCC +0.111, all better in 5 of 5 folds.

Full numbers are in `results/metrics/phase7/ddg_comparison.json` and `all_experiments_summary.csv`. Figures are in `results/figures/phase7/`.

## Interpretation

1. **ΔΔG carries real pathogenicity signal.** It is the top feature when only physicochemistry is available, and it raises ROC-AUC by about 0.07 there, with a CI well away from zero.
2. **That signal is almost entirely redundant with the stronger features.**
   - On the standard baseline, ΔΔG gives a gain of at most 0.002 ROC-AUC. Only for XGBoost does the CI clearly exclude zero. The gain is statistically detectable at best, and not practically meaningful.
   - Once AlphaFold structure features (pLDDT, burial, accessibility) are present, the gain is indistinguishable from zero for every model. ΔΔG mostly re-expresses "the residue is buried in a well-folded region".
3. **Without the circular predictors** (the intrinsic set), ΔΔG adds about +0.002 for random forest, consistent across all 5 folds but with a CI that touches zero.
4. **Performance did not drop** in any configuration. The ΔΔG feature is not harmful, only largely redundant.

Phase 8 examines why: ΔΔG distributions, the low-pLDDT regions, per-gene effects and confounders.
