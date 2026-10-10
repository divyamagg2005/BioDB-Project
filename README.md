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
├── src/biodb/          # pipeline package (config.py, common.py, one module per phase)
├── tests/              # pytest suite
├── data/
│   ├── raw/            # downloads and API caches (not versioned)
│   ├── processed/      # phaseN/ outputs (small, versioned)
│   └── final/
├── docs/               # per-phase reports
├── models/
├── results/
├── legacy/             # archived 5-gene pilot and first 100-gene attempt (see legacy/README.md)
├── pyproject.toml
└── requirements.txt
```

## Phases

Each phase ends with validation, tests, a reviewed diff, one commit, and a push.

| Phase | Objective | Status |
|---|---|---|
| 0 | Repository cleanup and shared configuration | Done |
| 1 | 100-gene cohort (25 pathogenic + 25 benign missense per gene) | Done |
| 2 | VEP annotation | Not started |
| 3 | AlphaFold structures and residue mapping | Not started |
| 4 | FoldX ΔΔG | Not started |
| 5 | Final ML dataset | Not started |
| 6 | Baseline model (conventional features) | Not started |
| 7 | Stability-enhanced model | Not started |
| 8 | Scientific comparison | Not started |
| 9 | Streamlit application | Not started |
| 10 | Deployment and documentation | Not started |

Shared settings (paths, random seed, assembly GRCh38, cohort size) live in `src/biodb/config.py`.

## Environment Setup

Create and activate a virtual environment from the project root:

### Windows PowerShell

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e .
```

### macOS/Linux

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install -e .
```

Verify the core imports:

```bash
python -c "import Bio, matplotlib, numpy, pandas, requests, scipy, seaborn, sklearn, xgboost; print('Core imports verified')"
```

Run the test suite with `python -m pytest`.
