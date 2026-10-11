# Phase 8: Results and Scientific Analysis

**Question:** Does mutation-induced protein stability information (FoldX ΔΔG on AlphaFold structures) improve the prediction of pathogenic versus benign missense variants?

**Short answer:** ΔΔG carries genuine pathogenicity signal, but it adds almost nothing to a model that already has standard variant-effect predictors or AlphaFold-derived structural context. Its information is largely the same "buried residue in a well-folded region" signal that those features already capture.

## Roles of each component

- **AlphaFold** predicts protein *structure*. It says nothing about pathogenicity directly. Here it supplied 3D models plus per-residue confidence (pLDDT), burial and secondary structure.
- **FoldX** estimates the *stability change* (ΔΔG, kcal/mol) a substitution causes in that predicted structure. It is a physics-based energy estimate, not a pathogenicity prediction.
- **The machine-learning models** predict *pathogenicity*, trained on ClinVar labels with gene-grouped cross-validation.

## Data and validation (recap)

- **Data:** 100 genes × (25 pathogenic + 25 benign) ClinVar missense variants (GRCh38, at least 1 review star, deduplicated), 5,000 in total.
- **Coverage:** every variant has VEP annotations, an AlphaFold residue and a FoldX ΔΔG.
- **Validation:** 5-fold cross-validation grouped by gene. Every gene is entirely in either training or test, so models are always scored on unseen proteins.
- **Statistics:** improvements are assessed with paired per-fold differences and a gene-level bootstrap 95% CI.

## Final comparison (random forest, ROC-AUC, mean over folds)

| Feature set | Without ΔΔG | With ΔΔG | Δ [95% CI] |
|---|---|---|---|
| Baseline: VEP predictors + gnomAD + substitution | 0.978 | 0.978 | +0.0008 [+0.0000, +0.0015] |
| Baseline + AlphaFold structure | 0.978 | 0.978 | +0.0004 [−0.0001, +0.0009] |
| Intrinsic: substitution + structure (no predictors, no gnomAD) | 0.886 | 0.888 | +0.0022 [−0.0001, +0.0045] |
| Substitution only | 0.742 | 0.809 | **+0.068 [+0.053, +0.083]** |

Logistic regression and XGBoost show the same pattern; full tables are in `docs/phase7_stability_report.md`. Baseline F1 is 0.925 and MCC 0.851; adding ΔΔG changes F1 by +0.001.

Figures (`results/figures/`):
- `phase6/roc_baseline.png`, `phase6/confusion_baseline.png`
- `phase7/roc_with_without_ddg.png`, `phase7/delta_auc_ci.png`
- `phase8/ddg_by_label_and_confidence.png`, `phase8/single_feature_auc.png`, `phase8/feature_importance_full_model.png`, `phase8/per_gene_ddg_auc.png`

## Why ΔΔG adds so little

### 1. ΔΔG is informative on its own

| Feature, alone | ROC-AUC |
|---|---|
| AlphaMissense | 0.952 |
| REVEL | 0.928 |
| CADD | 0.878 |
| PolyPhen-2 / SIFT | 0.854 / 0.839 |
| gnomAD absent | 0.835 |
| Relative solvent accessibility | 0.806 |
| Neighbour count (burial) | 0.802 |
| pLDDT | 0.753 |
| **FoldX ΔΔG** | **0.733** |
| GERP conservation | 0.732 |
| BLOSUM62 | 0.686 |

ΔΔG is comparable to conservation. However, it is weaker than the AlphaFold burial features computed from the same structures, and much weaker than the learned predictors.

### 2. ΔΔG works only where the structure is reliable

| Stratum | n | ΔΔG ROC-AUC | Median ΔΔG, pathogenic / benign |
|---|---:|---:|---|
| pLDDT ≥ 70 | 3,095 | 0.755 | 2.57 / 0.57 |
| pLDDT < 70 | 1,905 | 0.535 | 0.25 / 0.19 |
| pLDDT ≥ 70 and buried | 1,903 | 0.770 | — |
| pLDDT ≥ 70 and exposed | 1,192 | 0.613 | — |

In low-confidence (largely disordered) regions FoldX returns near-zero ΔΔG for both classes. 38% of the cohort, mostly benign variants (1,459 of 1,905), falls there.

Per gene, ΔΔG ROC-AUC has a median of 0.75. It exceeds 0.7 in 64 of 100 genes and is below 0.5 in 6. It rises with the gene's median pLDDT (Spearman 0.53).

### 3. ΔΔG overlaps with burial and the predictors

Spearman correlation of ΔΔG with:

| Feature | ρ |
|---|---:|
| relative solvent accessibility | −0.51 |
| neighbour count | 0.50 |
| AlphaMissense | 0.43 |
| pLDDT | 0.40 |
| REVEL | 0.39 |
| CADD | 0.32 |
| BLOSUM62 | −0.31 |

Destabilising mutations are mostly buried, well-folded residues. That's the same information RSA, neighbour count and pLDDT encode directly, and that AlphaMissense learned (AlphaMissense itself uses AlphaFold-derived structural context).

### 4. Where it helps a little

On the 900 variants where the baseline is uncertain (out-of-fold probability between 0.2 and 0.8), adding ΔΔG raises ROC-AUC from 0.769 to 0.787 and accuracy from 0.718 to 0.723. This subset is selected using the baseline's own predictions, so treat it as exploratory.

No gain is seen on the 293 variants where AlphaMissense itself is ambiguous (0.851 → 0.850).

In the full model, ΔΔG ranks 11th of 26 features by random-forest importance.

## Confounders and limitations

- **Circularity.** ClinVar curators use computational predictors (ACMG PP3/BP4) and population frequency (BA1/BS1) as evidence. The near-ceiling baseline (0.978) partly reflects labels that were derived from the same evidence. The "intrinsic" comparison removes these features; there ΔΔG still adds only about +0.002.
- **Cohort selection.** Genes were chosen for having at least 25 variants per class, 2,700 residues or fewer, and a full-length AlphaFold model, which favours well-studied disease genes. 64% of labels are 1-star (single submitter). Results may not carry over to other genes or stricter labels.
- **Location confound.** Benign variants are concentrated in disordered, exposed regions and pathogenic ones in folded cores. Every structural feature, ΔΔG included, partly detects *where* a variant is rather than *what* it does.
- **FoldX limitations.**
  - ΔΔG is a single static-structure estimate on predicted models.
  - Very large values (up to 45 kcal/mol) reflect steric clashes for buried small-to-large substitutions.
  - FoldX does not model loss of function through binding, catalysis or regulation; many pathogenic missense variants act that way, e.g. TP53 DNA-contact residues such as R248.
- **Fixed hyperparameters.** No tuning was done, so absolute numbers could move slightly with tuning. The paired comparisons use identical settings on both sides.
- **Class balance.** The cohort is balanced 50/50 by design. Real-world precision would be lower, because benign variants vastly outnumber pathogenic ones.

## Conclusion

In this 100-gene, 5,000-variant cohort, FoldX ΔΔG computed on AlphaFold structures is a meaningful but weak, structure-dependent signal of pathogenicity:
- ROC-AUC 0.73 alone, and 0.77 in confidently predicted buried residues;
- strong gains only when nothing better is available (+0.068 over physicochemistry alone).

When modern variant-effect predictors or simple AlphaFold burial and confidence features are present, ΔΔG is redundant: Δ ROC-AUC ≤ 0.002, not practically meaningful. Stability modelling may be most useful as interpretable mechanistic evidence for individual variants in well-folded domains, rather than as an accuracy booster for a classifier.

Run with `python -m biodb.phase8_analysis`. It reads saved outputs only, with no retraining. Numbers are in `results/metrics/phase8/analysis.json` and `per_gene_ddg_auc.csv`.
