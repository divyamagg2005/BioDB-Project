# Prediction of Pathogenic and Benign Genetic Variants Using Machine Learning

Class project by:

- Divyam - 23BCB0003
- Daksh Manchanda - 23BCB0047

> Research and teaching project. Not a clinical tool. Predictions must not be used for diagnosis or patient care.

---

## Contents

1. [What this project does](#1-what-this-project-does)
2. [Headline results](#2-headline-results)
3. [Quick start](#3-quick-start)
4. [The pipeline, phase by phase](#4-the-pipeline-phase-by-phase)
5. [Full results and scores](#5-full-results-and-scores)
6. [What the scores mean](#6-what-the-scores-mean)
7. [Worked examples](#7-worked-examples)
8. [The Variant Explorer app](#8-the-variant-explorer-app)
9. [Glossary of terms used in the app](#9-glossary-of-terms-used-in-the-app)
10. [Limitations](#10-limitations)
11. [Responsible use](#11-responsible-use)
12. [Reproducing everything](#12-reproducing-everything)
13. [Deploying the app](#13-deploying-the-app)
14. [Repository layout](#14-repository-layout)
15. [Future work](#15-future-work)

---

## 1. What this project does

A **missense variant** is a DNA change that swaps one amino acid in a protein for another (for example, tyrosine 220 to cysteine in p53, written `TP53 p.Y220C`). Some missense variants cause disease (*pathogenic*); most are harmless (*benign*). Clinical labs classify them using evidence such as population frequency, conservation and computer predictions.

This project builds a machine-learning pipeline that predicts **pathogenic vs benign** for missense variants from the ClinVar database. It asks one research question:

> **Does protein stability information, meaning how much a mutation destabilises the 3D protein structure, improve pathogenicity prediction beyond conventional variant features?**

To answer it, the project:

1. selects 5,000 high-quality ClinVar missense variants across 100 genes;
2. annotates them with standard predictors (VEP: CADD, REVEL, AlphaMissense and others);
3. maps each one onto its **AlphaFold** predicted protein structure;
4. computes the stability change **ΔΔG** for every variant with **FoldX**;
5. trains models with and without ΔΔG, using gene-held-out validation;
6. compares them statistically and explains the outcome;
7. presents everything in an interactive **Streamlit app**.

Three tools have three different jobs:

| Tool | What it does | What it does *not* do |
|---|---|---|
| **AlphaFold** | Predicts the 3D structure of a protein | Does not say whether a mutation is harmful |
| **FoldX** | Estimates how much a mutation changes protein stability (ΔΔG) | Does not predict disease |
| **Machine learning** | Predicts pathogenic vs benign from all features | Is only as good as its labels and features |

---

## 2. Headline results

| Result | Value |
|---|---|
| Best baseline model (random forest, VEP + substitution features) | **ROC-AUC 0.978**, accuracy 0.925, F1 0.925, MCC 0.851 |
| Same model with FoldX ΔΔG added | ROC-AUC 0.978, accuracy 0.926, F1 0.926, MCC 0.853 |
| Gain from ΔΔG on the baseline | **+0.0008 ROC-AUC** (95% CI +0.0000 to +0.0015). Real but negligible. |
| ΔΔG as the only feature | ROC-AUC 0.733 |
| Gain from ΔΔG when only amino-acid properties are available | **+0.068 ROC-AUC** (95% CI +0.053 to +0.083) |
| FoldX compute | 5,000/5,000 variants, 0 failures, 1 h 22 min on a 32-core AWS machine |

**Conclusion:** ΔΔG carries genuine pathogenicity signal, but that signal is almost entirely redundant with the features good models already have: modern predictors such as AlphaMissense, and simple AlphaFold measures of whether a residue is buried and confidently folded. Stability helps most where nothing better is available. Its main value is as an *explanation* for individual variants in well-folded protein regions.

---

## 3. Quick start

```bash
git clone https://github.com/divyamagg2005/BioDB-Project.git
cd BioDB-Project
python -m venv .venv
# Windows:      .\.venv\Scripts\Activate.ps1
# macOS/Linux:  source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .

python -m pytest                  # 70 tests
streamlit run streamlit_app.py    # opens the Variant Explorer at http://localhost:8501
```

Phase outputs, metrics, figures and the small models are versioned in the repository, so the app and Phases 5–9 work straight after cloning. You don't need FoldX, AWS or a network connection, except for loading 3D structures in the app.

---

## 4. The pipeline, phase by phase

Each phase is one git commit and one tag (`phase0` … `phase10`; Phase 1 is tagged `phase1-v2`) and has a detailed report in `docs/`.

### Phase 0: Cleanup and shared configuration

**Goal:** a clean, reproducible code base.

- Moved the original 5-gene pilot and a buggy first 100-gene attempt to `legacy/`, with reasons in [legacy/README.md](legacy/README.md).
- Created the `src/biodb/` package. `config.py` holds all paths, the random seed (42), the genome assembly (GRCh38) and the cohort size. `common.py` holds the shared parsers.
- Added `pyproject.toml` and a pytest suite.

### Phase 1: Building the 100-gene cohort

**Goal:** 100 genes × (25 pathogenic + 25 benign) = **5,000** reliable missense variants.

**Source:** ClinVar `variant_summary.txt.gz` (9,290,347 rows, downloaded 2026-10-10, checksum-verified).

The filters, in order:

| Step | Rule | Rows removed |
|---|---|---:|
| Assembly | GRCh38 only. ClinVar lists every variant once per assembly, so mixing assemblies duplicates rows. | 4,678,836 |
| Label | Pathogenic / Likely pathogenic → 1; Benign / Likely benign → 0. VUS and conflicting labels dropped. | 2,823,132 |
| Type | Single-nucleotide variant | 251,876 |
| Origin | Germline | 17,164 |
| Review status | At least 1 star ("criteria provided") | 44,303 |
| Missense | Exactly one amino-acid change between two standard amino acids | 1,276,031 |
| Start loss | `p.Met1` changes excluded (they stop translation rather than change a residue) | 1,664 |
| Gene | Unambiguous single gene | 829 |

Then:
- duplicates were removed: by VariationID, and by protein change within a gene;
- every gene was mapped to one reviewed UniProt protein with a full-length AlphaFold model (≤ 2,700 residues);
- every variant's reference amino acid was checked against that protein sequence.

157 genes qualified, 115 passed structure mapping, and the top 100 were taken. Within each gene, 25 variants per class were drawn at random (seed 42).

**Example:** `TP53 p.Y220C` (ClinVar VCV000127819) is in the cohort as pathogenic. Its tyrosine-220 matches residue 220 of UniProt P04637.

**Output:** `data/processed/phase1/cohort.csv`, 5,000 rows. All 18 validation checks pass, and a rerun reproduces the file byte for byte. Report: [docs/phase1_cohort_report.md](docs/phase1_cohort_report.md).

### Phase 2: VEP annotation

**Goal:** attach standard variant-effect predictors to every variant.

- Used the Ensembl **VEP** REST API, 25 batches of 200 variants each.
- For each variant, chose the transcript whose protein change matches ClinVar. MANE Select was preferred and was used for 4,998 of 5,000 variants.
- Extracted CADD, REVEL, AlphaMissense, SIFT, PolyPhen-2, GERP conservation and gnomAD allele frequency.
- The ClinVar fields VEP also returns are **never read**, so the label cannot leak into the features.

| Feature | Coverage | Median, pathogenic | Median, benign |
|---|---:|---:|---:|
| CADD | 100% | 27.7 | 21.6 |
| REVEL | 99.0% | 0.887 | 0.271 |
| AlphaMissense | 99.0% | 0.989 | 0.109 |
| SIFT | 100% | 0.00 | 0.06 |
| PolyPhen-2 | 97.0% | 0.990 | 0.038 |
| gnomAD AF present | 51.7% | 18% of pathogenic | 85% of benign |

**Example:** `TP53 p.Y220C` has AlphaMissense 0.852, REVEL 0.936, CADD 28.7, and gnomAD AF 6.6 × 10⁻⁶ (extremely rare).

Report: [docs/phase2_vep_report.md](docs/phase2_vep_report.md).

### Phase 3: AlphaFold structures

**Goal:** place every variant on its 3D protein structure and describe its surroundings.

- Downloaded 100 AlphaFold DB v6 models; each sequence is identical to UniProt.
- All 5,000 variants mapped to the correct residue.
- Computed, per residue: **pLDDT** (confidence), **relative solvent accessibility** (how exposed), **neighbour count** (how buried) and **secondary structure** (helix, strand or coil).

| | Pathogenic | Benign |
|---|---:|---:|
| Median pLDDT | 89.3 | 52.4 |
| Median relative solvent accessibility | 0.14 (buried) | 0.61 (exposed) |
| Median neighbour count | 19 | 8 |

Pathogenic variants sit mostly in buried, well-folded cores. Benign variants sit mostly in floppy, exposed regions.

**Example:** in TP53, R175 (a known structural hotspot) is buried (RSA 0.03). R248, a DNA-contact residue, is exposed (RSA 0.72).

Report: [docs/phase3_structure_report.md](docs/phase3_structure_report.md).

### Phase 4: FoldX ΔΔG (the heavy compute)

**Goal:** a stability change for every variant.

Steps:
1. Convert each AlphaFold model to PDB format.
2. Run FoldX `RepairPDB` once per protein.
3. Run FoldX `BuildModel` (3 runs averaged) once per variant.

The FoldX mutation format includes the chain, e.g. `YA220C`. ΔΔG is read from the `total energy` column and cross-checked against the wild-type and mutant energies.

- **Compute:** one AWS `c7a.8xlarge` spot instance (32 vCPU). 1 h 22 min of wall time, about 42 CPU-hours. 5,000/5,000 variants succeeded.
- **Sanity check:** `TP53 p.Y220C` gives ΔΔG **+3.48 kcal/mol**, close to the ~4 kcal/mol measured experimentally for this well-known destabilising mutation.

| Group | n | Median ΔΔG (kcal/mol) |
|---|---:|---:|
| Pathogenic, confident structure (pLDDT ≥ 70) | 2,054 | 2.57 |
| Benign, confident structure | 1,041 | 0.57 |
| Pathogenic, low confidence (pLDDT < 70) | 446 | 0.25 |
| Benign, low confidence | 1,459 | 0.19 |

Report: [docs/phase4_foldx_report.md](docs/phase4_foldx_report.md).

### Phase 5: Final ML dataset

**Goal:** one table, one row per variant, all features.

`data/processed/phase5/ml_dataset.csv` has 5,000 rows and 26 features in four groups:

| Group | Features |
|---|---|
| Conventional (8) | CADD, REVEL, AlphaMissense, SIFT, PolyPhen-2, GERP, gnomAD AF, gnomAD-absent flag |
| Substitution (11) | BLOSUM62, hydropathy change, volume change (signed and absolute), charge change, to/from proline, glycine and cysteine |
| Structure (6) | pLDDT, windowed pLDDT, relative solvent accessibility, neighbour count, helix, strand |
| Stability (1) | FoldX ΔΔG |

Missing values stay missing in the file; they are filled in inside model training, using training data only. No ClinVar-derived field is a feature.

Report: [docs/phase5_dataset_report.md](docs/phase5_dataset_report.md).

### Phase 6: Baseline model

**Goal:** the control experiment, using conventional and substitution features only (19 features, no structure, no ΔΔG).

- **Validation:** 5-fold cross-validation **grouped by gene**. Each fold tests on 20 genes the model never saw during training.
- **Models:** logistic regression, random forest and XGBoost.
- **Result:** random forest is best at ROC-AUC **0.978 ± 0.005**. See section 5 for every metric.
- **Most important features:** AlphaMissense, REVEL and gnomAD absence.

Report: [docs/phase6_baseline_report.md](docs/phase6_baseline_report.md).

### Phase 7: Stability-enhanced models

**Goal:** add ΔΔG and measure the difference fairly. The data, folds, models and settings are identical to Phase 6; only the feature set changes.

ΔΔG is added to four reference feature sets:

| Reference | Question | ROC-AUC change from ΔΔG, random forest [95% CI] |
|---|---|---|
| Baseline | Does ΔΔG improve the standard model? | +0.0008 [+0.0000, +0.0015] |
| Baseline + structure | Does ΔΔG add anything beyond burial and confidence? | +0.0004 [−0.0001, +0.0009] |
| Intrinsic (substitution + structure, no predictors) | Does it help without the predictors? | +0.0022 [−0.0001, +0.0045] |
| Substitution only | Does ΔΔG carry signal at all? | **+0.068 [+0.053, +0.083]** |

The confidence intervals come from a gene-level bootstrap: genes are resampled 1,000 times.

Report: [docs/phase7_stability_report.md](docs/phase7_stability_report.md).

### Phase 8: Scientific analysis

**Goal:** explain *why* ΔΔG adds so little.

- **It works only on reliable structure.** ΔΔG ROC-AUC is 0.755 where pLDDT ≥ 70, 0.770 in confident buried residues, and 0.535 (near random) where pLDDT < 70.
- **It overlaps with burial.** Spearman correlation with relative solvent accessibility is −0.51 and with neighbour count +0.50. Destabilising mutations are the buried ones, which those features already capture.
- **It overlaps with AlphaMissense** (ρ = 0.43). AlphaMissense itself learned from AlphaFold-derived structural context.
- **It is consistent across genes.** The median per-gene ΔΔG ROC-AUC is 0.75, and 64 of 100 genes are above 0.7. The per-gene value rises with structure confidence.
- **It helps a little on hard cases.** On the 900 variants where the baseline is uncertain, adding ΔΔG raises ROC-AUC from 0.769 to 0.787.

Report: [docs/phase8_results.md](docs/phase8_results.md).

### Phase 9: Streamlit Variant Explorer

**Goal:** an interactive app showing predictions, ΔΔG and the 3D structure for any cohort variant. See [section 8](#8-the-variant-explorer-app). Report: [docs/phase9_app.md](docs/phase9_app.md).

### Phase 10: Release

**Goal:** a project anyone can install and verify.

- **Secret scan:** no keys, tokens or FoldX binary anywhere in the files or the git history.
- **Size:** the largest tracked file is 1.3 MB.
- **Library versions:** scikit-learn and XGBoost are pinned to the versions the saved models were trained with.
- **Clean install:** a fresh virtual environment passes all 70 tests.

---

## 5. Full results and scores

All numbers are means over the 5 gene-grouped cross-validation folds. Each fold is tested on 20 genes the model never saw. The threshold for accuracy, precision, recall, F1 and MCC is 0.5.

| Feature set | Model | Accuracy | Precision | Recall | F1 | ROC-AUC | PR-AUC | MCC | Brier |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **Baseline** | Logistic regression | 0.913 | 0.909 | 0.918 | 0.914 | 0.974 | 0.973 | 0.826 | 0.062 |
| | Random forest | **0.925** | 0.928 | 0.922 | **0.925** | **0.978** | 0.974 | **0.851** | 0.058 |
| | XGBoost | 0.919 | 0.915 | 0.925 | 0.920 | 0.974 | 0.972 | 0.839 | 0.063 |
| **Baseline + ΔΔG** | Logistic regression | 0.915 | 0.913 | 0.919 | 0.916 | 0.975 | 0.974 | 0.831 | 0.061 |
| | Random forest | 0.926 | 0.929 | 0.924 | 0.926 | 0.978 | 0.976 | 0.853 | 0.056 |
| | XGBoost | 0.921 | 0.917 | 0.927 | 0.922 | 0.976 | 0.975 | 0.843 | 0.060 |
| **Baseline + structure** | Logistic regression | 0.917 | 0.912 | 0.924 | 0.918 | 0.975 | 0.974 | 0.834 | 0.061 |
| | Random forest | 0.925 | 0.928 | 0.922 | 0.924 | 0.978 | 0.975 | 0.850 | 0.056 |
| | XGBoost | 0.923 | 0.919 | 0.930 | 0.924 | 0.976 | 0.976 | 0.847 | 0.059 |
| **Baseline + structure + ΔΔG** | Logistic regression | 0.917 | 0.913 | 0.922 | 0.917 | 0.975 | 0.974 | 0.834 | 0.061 |
| | Random forest | **0.928** | **0.931** | 0.925 | **0.928** | 0.978 | 0.975 | **0.856** | 0.056 |
| | XGBoost | 0.923 | 0.919 | 0.928 | 0.924 | 0.977 | 0.976 | 0.847 | 0.058 |
| **Intrinsic** (substitution + structure) | Random forest | 0.809 | 0.802 | 0.823 | 0.812 | 0.886 | 0.882 | 0.619 | 0.137 |
| **Intrinsic + ΔΔG** | Random forest | 0.808 | 0.804 | 0.818 | 0.810 | 0.888 | 0.888 | 0.618 | 0.136 |
| **Substitution only** | Random forest | 0.675 | 0.689 | 0.638 | 0.662 | 0.742 | 0.745 | 0.350 | 0.206 |
| **Substitution + ΔΔG** | Random forest | 0.730 | 0.748 | 0.694 | 0.719 | 0.809 | 0.826 | 0.461 | 0.178 |
| | XGBoost | 0.740 | 0.767 | 0.690 | 0.726 | 0.816 | 0.833 | 0.483 | 0.175 |

Baseline random forest confusion matrix (all 5,000 variants, out-of-fold):

| | Predicted benign | Predicted pathogenic |
|---|---:|---:|
| **Actually benign** | 2,320 (true negative) | 180 (false positive) |
| **Actually pathogenic** | 194 (false negative) | 2,306 (true positive) |

How well each feature separates the classes on its own (ROC-AUC): AlphaMissense 0.952, REVEL 0.928, CADD 0.878, PolyPhen-2 0.854, SIFT 0.839, gnomAD-absent 0.835, solvent accessibility 0.806, neighbour count 0.802, pLDDT 0.753, **ΔΔG 0.733**, GERP 0.732, BLOSUM62 0.686.

Every number is in `results/metrics/` and every figure in `results/figures/`.

---

## 6. What the scores mean

All metrics compare the model's prediction with the ClinVar label. "Positive" means pathogenic.

| Metric | Plain meaning | Range (best) |
|---|---|---|
| **Accuracy** | Fraction of variants classified correctly | 0–1 (1) |
| **Precision** | Of the variants called pathogenic, the fraction that really are | 0–1 (1) |
| **Recall** (sensitivity) | Of the truly pathogenic variants, the fraction caught | 0–1 (1) |
| **F1** | Balance of precision and recall (their harmonic mean) | 0–1 (1) |
| **ROC-AUC** | Probability that a random pathogenic variant gets a higher score than a random benign one. 0.5 is a coin flip. It doesn't depend on the 0.5 threshold. | 0.5–1 (1) |
| **PR-AUC** | Area under the precision–recall curve; summarises precision across all thresholds | 0–1 (1) |
| **MCC** (Matthews correlation) | Correlation between predictions and truth that uses all four confusion-matrix cells. 0 is random. | −1 to 1 (1) |
| **Brier score** | Mean squared error of the predicted probabilities. Measures calibration: are "80%" predictions right 80% of the time? | 0–1 (**0**) |
| **± value** | Standard deviation across the 5 folds, i.e. how stable the score is from one group of genes to another | smaller is steadier |
| **95% CI** | Range that plausibly contains the true improvement. If it excludes 0, the improvement is unlikely to be noise. | — |

Other evaluation terms:
- **Gene-grouped cross-validation:** the 100 genes are split into 5 groups. The model trains on 4 groups and is tested on the 5th, which it never saw, and this repeats 5 times. This stops the model from "memorising" genes and shows how it would do on a new gene.
- **Out-of-fold prediction:** the prediction each variant got while its gene was in the test group. Every variant has exactly one.
- **Gene-level bootstrap:** the whole comparison is recomputed 1,000 times on random resamples of the genes. The spread of results gives the 95% CI.

---

## 7. Worked examples

All probabilities below are out-of-fold XGBoost predictions: the model never saw any variant of that gene.

### Example 1: TP53 p.Y220C (pathogenic, correctly predicted)

| Feature | Value | Interpretation |
|---|---|---|
| AlphaMissense / REVEL / CADD | 0.852 / 0.936 / 28.7 | All predict damaging |
| gnomAD AF | 6.6 × 10⁻⁶ | Extremely rare in the population |
| pLDDT | 97.7 | Very confident structure |
| Relative solvent accessibility / neighbours | 0.10 / 23 | Buried in the protein core |
| FoldX ΔΔG | **+3.48 kcal/mol** | Strongly destabilising (measured: ~4) |
| **Baseline → with ΔΔG** | **0.814 → 0.860** pathogenic | Both predict pathogenic; ΔΔG adds confidence |

This is a textbook "destabilising cavity" mutation: a buried tyrosine replaced by the much smaller cysteine.

### Example 2: TP53 p.P8L (benign, correctly predicted)

| Feature | Value | Interpretation |
|---|---|---|
| AlphaMissense / REVEL / CADD | 0.055 / 0.514 / 12.3 | Mostly predict tolerated |
| pLDDT | 39.3 | Disordered tail, no stable structure |
| Relative solvent accessibility / neighbours | 0.71 / 4 | Exposed, almost no contacts |
| FoldX ΔΔG | −0.37 | Near zero; not meaningful at this confidence |
| **Baseline → with ΔΔG** | **0.062 → 0.047** pathogenic | Both predict benign |

### Example 3: BRCA1 p.A1708E (where ΔΔG changes the answer)

| Feature | Value | Interpretation |
|---|---|---|
| AlphaMissense / REVEL / CADD | 0.971 / 0.831 / 25.5 | Damaging, but the baseline was still unsure |
| pLDDT / relative solvent accessibility | 96.1 / 0.02 | Confident, completely buried (BRCT domain) |
| FoldX ΔΔG | **+12.0 kcal/mol** | Severe destabilisation: a charged glutamate forced into the core |
| **Baseline → with ΔΔG** | **0.283 → 0.780** pathogenic | **Wrong (benign) → right (pathogenic)** |

BRCA1 A1708E is a well-studied BRCT-domain variant known to misfold the domain. This is exactly the kind of case where stability is informative. Such cases are rare, which is why the average gain is small.

### Example 4: ABCA12 p.Q515H (a limitation)

The label is pathogenic, but AlphaMissense (0.095) and REVEL (0.136) both call it benign. The residue sits in a low-confidence region (pLDDT 50), and ΔΔG is small (0.83). Both models predict benign (0.003). Without strong evidence in any feature, no model can recover such variants; they may act through mechanisms (splicing, binding, regulation) that none of these features capture.

---

## 8. The Variant Explorer app

```bash
streamlit run streamlit_app.py
```

**How to use it:**
1. Choose a **Gene** in the sidebar (default TP53).
2. Choose a **Variant** from the list, or type a protein change such as `Y220C` in the box.
3. Read the page top to bottom.

**What the page shows:**

| Section | Content |
|---|---|
| Header | Gene, protein change, and links to UniProt, AlphaFold DB and ClinVar |
| Machine learning: pathogenicity prediction | Predicted class and pathogenic/benign probability from two models: **Baseline (VEP + substitution)** and **Baseline + FoldX ΔΔG** |
| Final model (expander) | Probabilities from models retrained on all 5,000 variants. Optimistic, because this variant was in their training data. |
| ClinVar classification (expander) | The real ClinVar label and review status: the "answer key" the model learned from |
| AlphaFold: predicted structure | Interactive 3D cartoon: drag to rotate, scroll *over the viewer* to zoom. The mutated residue is shown in magenta. |
| FoldX: estimated stability change | ΔΔG in kcal/mol, with a warning when the structure is low-confidence |
| Structural context | pLDDT, relative solvent accessibility, neighbour count, secondary structure |
| Conventional features (VEP) | AlphaMissense, REVEL, CADD, SIFT, PolyPhen-2, GERP, gnomAD AF, BLOSUM62 |

Only the 5,000 cohort variants are supported, because each needs precomputed VEP and FoldX results. Any other change shows an error message. The app never retrains models.

---

## 9. Glossary of terms used in the app

### Variant and database terms

| Term | Meaning |
|---|---|
| **p.Y220C** | Protein change notation: amino acid **Y** (tyrosine) at position **220** becomes **C** (cysteine). One-letter amino-acid codes are used. |
| **Missense variant** | A DNA change that substitutes one amino acid for another |
| **ClinVar** | NCBI's public database of variant classifications submitted by clinical labs |
| **Pathogenic / Benign** | Disease-causing / harmless. "Likely" versions are included in each class. |
| **VCV number** | ClinVar's variation ID, e.g. VCV000127819. It links to the ClinVar record. |
| **Review status (stars)** | How well supported the ClinVar classification is: 1 star = single submitter with criteria; 2 = multiple agreeing submitters; 3 = expert panel; 4 = practice guideline |
| **UniProt** | The reference protein database. Its accession (e.g. P04637 for p53) identifies the exact protein sequence. |
| **GRCh38** | The current human reference genome assembly used for all coordinates |

### Prediction terms

| Term | Meaning |
|---|---|
| **Pathogenic probability** | The model's estimated chance that the variant is pathogenic (0–1). Benign probability = 1 − pathogenic probability. |
| **Predicted class** | Pathogenic if the pathogenic probability is ≥ 0.5, otherwise Benign |
| **Out-of-fold** | The probability came from a model trained *without any variant of this gene*. This is the honest estimate shown by default. |
| **Final model** | The model retrained on all 5,000 variants. Its numbers for cohort variants are optimistic. |
| **Baseline model** | Uses only conventional (VEP) and substitution features |
| **Baseline + FoldX ΔΔG** | The same model with the stability feature added |
| **XGBoost** | Gradient-boosted decision trees, the model type used in the app |

### AlphaFold / structure terms

| Term | Meaning |
|---|---|
| **AlphaFold** | DeepMind/EMBL-EBI AI system that predicts protein 3D structures. AlphaFold DB stores predictions for most human proteins. |
| **AF-P04637-F1** | AlphaFold DB entry ID: UniProt accession P04637, fragment 1 (the full protein) |
| **pLDDT** | AlphaFold's per-residue confidence (0–100). **> 90** very high (dark blue), **70–90** confident (light blue), **50–70** low (yellow), **< 50** very low, usually disordered (orange). |
| **Relative solvent accessibility (RSA)** | How exposed the residue is to the surrounding water: 0 = fully buried, about 1 = fully exposed. Buried residues are usually more sensitive to mutation. |
| **Neighbour count (10 Å)** | Number of residues within 10 ångström (1 Å = 0.1 nm) of this residue. High = packed inside the protein; low = on the surface or in a floppy region. |
| **Secondary structure** | Local backbone shape: **helix** (α-helix), **strand** (β-sheet) or **coil** (loops and turns) |
| **Magenta residue** | The mutated amino acid, drawn as sticks with a sphere on its central carbon (Cα) |

### FoldX / stability terms

| Term | Meaning |
|---|---|
| **FoldX** | Software that estimates protein stability from a 3D structure using an empirical energy function |
| **ΔΔG (kcal/mol)** | Change in folding free energy caused by the mutation: ΔG(mutant) − ΔG(wild type). **Positive = destabilising.** Roughly: < 1 neutral, 1–2 mildly destabilising, > 2 strongly destabilising. |
| **Low-confidence warning** | Shown when pLDDT < 70: the structure there is unreliable, so ΔΔG is too |
| **RepairPDB / BuildModel** | FoldX commands: first optimise the structure, then build the mutant and compute energies |

### VEP feature terms

| Term | Meaning | Direction |
|---|---|---|
| **VEP** | Ensembl Variant Effect Predictor: annotates variants with consequences and scores | — |
| **AlphaMissense** | DeepMind model predicting missense pathogenicity (0–1); < 0.34 likely benign, > 0.564 likely pathogenic | higher = more damaging |
| **REVEL** | Ensemble of 13 predictors (0–1); > 0.5 suggests pathogenic | higher = more damaging |
| **CADD (phred)** | Genome-wide deleteriousness score; 20 = top 1%, 30 = top 0.1% of possible variants | higher = more damaging |
| **SIFT** | Predicts whether a substitution is tolerated (0–1); < 0.05 = deleterious | **lower** = more damaging |
| **PolyPhen-2** | Predicts impact on protein function (0–1); ≥ 0.909 = probably damaging (HumVar model used by VEP) | higher = more damaging |
| **GERP** | Evolutionary conservation of the DNA position; > 2 = conserved | higher = more conserved |
| **gnomAD AF** | Allele frequency in ~800,000 people (gnomAD). "Absent from gnomAD" means never observed; disease variants are usually very rare. | lower / absent = more suspicious |
| **BLOSUM62** | How often this amino-acid swap is seen in evolution; negative = unusual swap | lower = more radical change |

---

## 10. Limitations

- **Circular evidence.** ClinVar curators use computational predictors (ACMG criteria PP3/BP4) and population frequency (BA1/BS1) when classifying variants. The near-perfect baseline is therefore partly re-learning the curators' own evidence.
- **Cohort bias.** Genes were selected for having at least 25 variants in each class and a full-length AlphaFold model, which favours well-studied disease genes. 64% of labels are 1-star (single submitter).
- **Location confound.** Benign variants cluster in disordered, exposed regions and pathogenic ones in folded cores. Structural features (ΔΔG included) partly detect *where* a variant is rather than *what* it does.
- **FoldX scope.** ΔΔG is estimated from one static, predicted structure of one protein chain. It misses mechanisms such as binding, catalysis, splicing and regulation, and is unreliable where pLDDT < 70.
- **Balanced classes.** The 50/50 design overstates real-world precision; in reality benign variants vastly outnumber pathogenic ones.
- **Fixed hyperparameters.** Models were not tuned. Comparisons are fair because both sides use identical settings, but absolute numbers could shift slightly with tuning.

---

## 11. Responsible use

- **Research and teaching only.** Not validated for clinical use.
- **Public data only:** ClinVar, Ensembl VEP, gnomAD summary frequencies, UniProt, AlphaFold DB. No patient data.
- **Honest probabilities in the app.** It shows out-of-fold probabilities by default; in-sample numbers are clearly labelled.
- **Inherited bias.** The models inherit ClinVar's biases toward well-studied genes and certain populations. Expect lower performance elsewhere.

---

## 12. Reproducing everything

```bash
python -m biodb.phase1_cohort          # needs data/raw/clinvar/variant_summary.txt.gz (NCBI FTP)
python -m biodb.phase2_vep             # Ensembl VEP REST, ~1 hour
python -m biodb.phase3_structure       # downloads 100 AlphaFold models (~126 MB)
python -m biodb.phase4_foldx all       # needs a licensed FoldX 5.1 binary (FOLDX_BIN); ~42 CPU-hours
python -m biodb.phase5_dataset
python -m biodb.phase6_baseline
python -m biodb.phase7_stability
python -m biodb.phase8_analysis
```

- **Steps 1–4** need network access and save raw downloads to `data/raw/` (not versioned). Steps 1 and 2 can rerun offline from their caches (`--offline`).
- **Step 4** ran on AWS EC2 (`c7a.8xlarge` spot, Ubuntu 24.04). FoldX is licensed separately (academic licence) and is not included. The FoldX 5.1 build used stops working on 31 Dec 2026.
- **Steps 5–8** take minutes and need only the versioned outputs.

---

## 13. Deploying the app

The app needs no secrets and no cloud compute.

1. Sign in at <https://share.streamlit.io> with the GitHub account that owns this repository.
2. **Create app**, then choose repository `divyamagg2005/BioDB-Project`, branch `main`, file `streamlit_app.py`.
3. Under **Advanced settings**, choose Python 3.11.
4. Click **Deploy**.

---

## 14. Repository layout

```text
src/biodb/             pipeline package
  config.py            paths, seed, assembly, cohort size
  common.py            shared parsers and helpers
  phase1_cohort.py … phase8_analysis.py   one module per phase
  modeling.py          shared training/evaluation (folds, models, metrics)
  app_data.py          data access for the app
streamlit_app.py       Variant Explorer app
tests/                 70 pytest tests and real FoldX output fixtures
data/processed/        versioned outputs of phases 1–5
results/metrics/       metrics, folds, out-of-fold predictions (phases 6–8)
results/figures/       ROC curves, confusion matrix, ΔΔG plots
models/                saved logistic-regression and XGBoost models
docs/                  detailed report per phase
legacy/                archived earlier attempts
```

---

## 15. Future work

- Add more genes and stricter labels (2+ review stars) to test robustness.
- Add mechanism-specific features (interface and binding-site proximity, catalytic residues) that monomer stability cannot capture.
- Tune hyperparameters with nested gene-grouped cross-validation, and calibrate probabilities to realistic class frequencies.
- Validate on variants newly classified in a later ClinVar release.
