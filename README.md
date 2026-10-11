# Prediction of Pathogenic and Benign Genetic Variants Using Machine Learning

Class project by:

- Divyam - 23BCB0003
- Daksh Manchanda - 23BCB0047

## 1. Problem

Most missense variants found in patients are of uncertain significance. This project builds an interpretable machine-learning pipeline that classifies ClinVar missense variants as pathogenic or benign. It asks one question:

> Does mutation-induced protein stability, computed with FoldX on AlphaFold structures, add useful signal beyond conventional variant-level features?

This is a class project, not a clinical tool. It deliberately avoids deep learning, molecular dynamics and genome-scale downloads.

**Answer, in short:** ΔΔG carries real signal (ROC-AUC 0.73 on its own, +0.068 over physicochemistry alone). However, it is almost fully redundant once standard predictors or AlphaFold burial and confidence features are present (Δ ROC-AUC ≤ 0.002). See [docs/phase8_results.md](docs/phase8_results.md).

## 2. Pipeline

| Phase | What | Output | Report |
|---|---|---|---|
| 0 | Repository cleanup, shared config | `src/biodb/config.py` | [legacy/README.md](legacy/README.md) |
| 1 | 100-gene ClinVar cohort, 25 pathogenic + 25 benign missense per gene | `data/processed/phase1/` | [phase1](docs/phase1_cohort_report.md) |
| 2 | VEP annotation (CADD, REVEL, AlphaMissense, SIFT, PolyPhen-2, GERP, gnomAD) | `data/processed/phase2/` | [phase2](docs/phase2_vep_report.md) |
| 3 | AlphaFold structures, residue mapping, structural features | `data/processed/phase3/` | [phase3](docs/phase3_structure_report.md) |
| 4 | FoldX ΔΔG for all 5,000 variants (AWS EC2) | `data/processed/phase4/` | [phase4](docs/phase4_foldx_report.md) |
| 5 | Final ML dataset, 26 features | `data/processed/phase5/` | [phase5](docs/phase5_dataset_report.md) |
| 6 | Baseline model (conventional + substitution features) | `results/*/phase6/`, `models/phase6/` | [phase6](docs/phase6_baseline_report.md) |
| 7 | Stability-enhanced models (+ ΔΔG), paired comparisons | `results/*/phase7/`, `models/phase7/` | [phase7](docs/phase7_stability_report.md) |
| 8 | Scientific analysis | `results/*/phase8/` | [phase8](docs/phase8_results.md) |
| 9 | Streamlit variant explorer | `streamlit_app.py` | [phase9](docs/phase9_app.md) |
| 10 | Release checks, documentation, deployment | this README | — |

Each phase is one commit and one git tag (`phase0` … `phase10`; Phase 1 is tagged `phase1-v2`).

## 3. Data

- **Labels:** ClinVar `variant_summary.txt.gz` (downloaded 2026-10-10; SHA-256 recorded in `data/processed/phase1/metadata.json`).
  - Pathogenic, Likely pathogenic and Pathogenic/Likely pathogenic map to `1`.
  - Benign, Likely benign and Benign/Likely benign map to `0`.
  - VUS, conflicting and mixed labels are excluded.
- **Cohort** (Phase 1):
  - variants: GRCh38, germline, single-nucleotide missense, at least 1 review star, deduplicated by VariationID and by protein change, start-loss changes excluded;
  - genes: mapped to one reviewed UniProt entry and a full-length AlphaFold model (≤ 2,700 residues), with the reference residue checked;
  - size: 100 genes, 50 variants each (25 per class), sampled at random with seed 42, giving 5,000 variants.
- **Annotations:** Ensembl VEP REST, using the MANE Select transcript. ClinVar fields returned by VEP are never read.
- **Structures:** AlphaFold DB v6 canonical models. They are downloaded on demand and not stored in git; checksums are in `data/processed/phase3/structures.csv`.
- **Stability:** FoldX 5.1. RepairPDB runs once per protein, then BuildModel with 3 runs per variant. ΔΔG is mutant minus wild type, in kcal/mol. The FoldX binary and licence are not in this repository.

## 4. Features (Phase 5)

| Group | Features |
|---|---|
| Conventional | CADD, REVEL, AlphaMissense, SIFT, PolyPhen-2, GERP, gnomAD AF, gnomAD-absent flag |
| Substitution | BLOSUM62, hydropathy change, volume change (signed and absolute), charge change, Pro/Gly/Cys flags |
| Structure (AlphaFold) | pLDDT, windowed pLDDT, relative solvent accessibility, CB neighbour count, helix/strand |
| Stability (FoldX) | ΔΔG |

Missing values stay missing in the dataset. Imputation is fitted inside each training fold.

## 5. Validation

- **Folds:** 5-fold `StratifiedGroupKFold` **grouped by gene** (20 genes per fold), saved once in `results/metrics/folds.csv` and reused by every experiment. Every model is always scored on genes it never saw.
- **Models:** logistic regression, random forest and XGBoost, with fixed hyperparameters and identical settings across experiments.
- **Metrics:** accuracy, precision, recall, F1, ROC-AUC, PR-AUC, MCC and Brier score (mean ± SD over folds).
- **Comparisons:** ΔΔG comparisons are paired. A gene-level bootstrap gives a 95% CI for the ROC-AUC difference.

## 6. Results

| Feature set (random forest) | ROC-AUC without ΔΔG | with ΔΔG | Δ [95% CI] |
|---|---|---|---|
| Baseline (VEP + gnomAD + substitution) | 0.978 | 0.978 | +0.0008 [+0.0000, +0.0015] |
| Baseline + AlphaFold structure | 0.978 | 0.978 | +0.0004 [−0.0001, +0.0009] |
| Intrinsic (substitution + structure) | 0.886 | 0.888 | +0.0022 [−0.0001, +0.0045] |
| Substitution only | 0.742 | 0.809 | +0.068 [+0.053, +0.083] |

- The baseline random forest reaches F1 0.925 and MCC 0.851.
- ΔΔG alone has ROC-AUC 0.733: 0.755 where pLDDT ≥ 70, and 0.535 where pLDDT < 70.
- ΔΔG correlates with burial (relative solvent accessibility ρ = −0.51, neighbour count ρ = 0.50) and with AlphaMissense (ρ = 0.43).

## 7. Limitations

- **Circularity.** ClinVar curators use computational predictors and population frequency as evidence, which inflates the baseline.
- **Cohort bias.** Gene selection favours well-studied disease genes. 64% of labels are 1-star.
- **Location confound.** Benign variants are concentrated in disordered, exposed regions.
- **FoldX scope.** FoldX gives a static-structure estimate. It does not capture binding, catalytic or regulatory mechanisms. Its values are unreliable where pLDDT < 70.
- **Class balance.** The 50/50 design overstates real-world precision.

Details are in [docs/phase8_results.md](docs/phase8_results.md).

## 8. Responsible use

- **Research and teaching only.** Predictions must not be used for diagnosis, reporting or clinical decisions.
- **Public data only.** All data are public: ClinVar, Ensembl VEP, gnomAD summaries, UniProt and AlphaFold DB. No patient-level or identifiable data are used.
- **Honest probabilities in the app.** The app shows out-of-fold probabilities, from models that never saw that gene. In-sample probabilities from the final models appear only in a clearly labelled section.
- **Known blind spots.** The models inherit ClinVar's ascertainment biases (populations, well-studied genes, submitter practices), and performance on other genes may be lower.

## 9. Setup and reproduction

```bash
python -m venv .venv
# Windows: .\.venv\Scripts\Activate.ps1    macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip install -e .
python -m pytest            # 70 tests: unit tests plus checks on the versioned outputs
```

Rebuilding from scratch needs network access. Phase 4 also needs a licensed FoldX 5.1 binary (`FOLDX_BIN`).

```bash
python -m biodb.phase1_cohort          # needs data/raw/clinvar/variant_summary.txt.gz
python -m biodb.phase2_vep
python -m biodb.phase3_structure
python -m biodb.phase4_foldx all       # ~42 CPU-hours; run on a many-core Linux machine
python -m biodb.phase5_dataset
python -m biodb.phase6_baseline
python -m biodb.phase7_stability
python -m biodb.phase8_analysis
```

Phases 5–9 run in minutes from the versioned outputs, without FoldX or network access.

## 10. Variant explorer app

```bash
streamlit run streamlit_app.py
```

Pick a gene and a cohort variant, or type a change such as `Y220C`. The app shows:
- baseline and baseline + ΔΔG pathogenicity probabilities;
- FoldX ΔΔG;
- the AlphaFold 3D structure, coloured by pLDDT, with the mutated residue highlighted;
- the VEP features.

### Deploying to Streamlit Community Cloud

The app needs no secrets and no AWS. It reads versioned files and fetches structures from AlphaFold DB.

1. Sign in at <https://share.streamlit.io> with the GitHub account that owns `divyamagg2005/BioDB-Project`.
2. **Create app**, then choose this repository, branch `main`, main file `streamlit_app.py`.
3. Under **Advanced settings**, choose Python 3.11 to match development.
4. Click **Deploy**. Dependencies install from `requirements.txt`.

## 11. Future work

- Add more genes and stricter labels (2+ review stars) to test robustness.
- Add mechanism-specific features (interface and binding-site proximity, catalytic residues) that FoldX monomer stability cannot capture.
- Tune hyperparameters with nested gene-grouped cross-validation, and calibrate probabilities for realistic class priors.
- Evaluate on an independent set, such as variants newly classified in a later ClinVar release.

## Repository layout

```text
src/biodb/          pipeline package: one module per phase, shared config/modeling
tests/              pytest suite (unit tests + output checks) and FoldX fixtures
data/processed/     versioned phase outputs (data/raw is not versioned)
results/            metrics, out-of-fold predictions, figures
models/             saved logistic-regression and XGBoost models
docs/               per-phase reports
streamlit_app.py    variant explorer
legacy/             archived 5-gene pilot and first 100-gene attempt
```
