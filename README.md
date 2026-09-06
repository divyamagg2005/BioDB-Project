# Prediction of Pathogenic and Benign Genetic Variants Using Machine Learning

Class project by:

- Divyam - 23BCB0003
- Daksh Manchanda - 23BCB0047

## Project Scope

This project will develop an interpretable machine learning pipeline for predicting whether a genetic variant is pathogenic or benign. The main investigation is whether protein structural and stability information derived from AlphaFold structures and mutation analysis adds useful signal beyond conventional variant-level features.

This is a class project, not a research-grade clinical system. The implementation will remain computationally manageable and will avoid unnecessary deep learning, molecular dynamics, transformers, and large database downloads.

## Data Policy

- ClinVar is the primary source for clinical variant labels.
- AlphaFold DB is the source for protein structures.
- Additional annotation sources may be used only when scientifically justified and documented.
- The project will not download the entire AlphaFold DB or dbNSFP database.
- Random Kaggle datasets will not be used as substitutes for ClinVar.

### Target Definition

- Pathogenic and Likely pathogenic: `1`
- Benign and Likely benign: `0`
- VUS and conflicting classifications: excluded

## Repository Structure

```text
.
├── data/
│   ├── raw/
│   ├── processed/
│   └── final/
├── docs/
├── models/
├── notebooks/
├── results/
│   ├── figures/
│   └── metrics/
├── src/
├── tests/
├── .gitignore
├── README.md
└── requirements.txt
```

## Phase 0: Environment Setup

Create and activate a virtual environment from the project root:

### Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

### macOS/Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Verify the core imports:

```bash
python -c "import Bio, matplotlib, numpy, pandas, requests, scipy, seaborn, sklearn, xgboost; print('Core imports verified')"
```

Phase 0 intentionally does not implement the machine learning pipeline. Data acquisition, preprocessing, feature engineering, model training, evaluation, and structural analysis will be added in later phases.
