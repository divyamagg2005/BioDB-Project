# Phase 9: Streamlit Variant Explorer

Run with `streamlit run streamlit_app.py` from the repository root. The data access code is in `src/biodb/app_data.py`, and the tests in `tests/test_app.py`.

## What it shows

For a selected cohort variant (gene drop-down, then a variant or a typed change such as `Y220C`):

| Section | Source | Meaning |
|---|---|---|
| Machine learning | `results/metrics/phase7/*_oof_predictions.csv` (XGBoost) | Pathogenic and benign probability, plus the predicted class at threshold 0.5, from the baseline model and the baseline + ΔΔG model |
| AlphaFold | AlphaFold DB mmCIF fetched at view time, rendered with 3Dmol.js | Predicted structure coloured by pLDDT, mutated residue in magenta |
| FoldX | `data/processed/phase5/ml_dataset.csv` | ΔΔG in kcal/mol, with a warning when pLDDT < 70 |
| Structural context | Phase 3 features | pLDDT, relative solvent accessibility, neighbour count, secondary structure |
| Conventional features | Phase 2 VEP annotations | AlphaMissense, REVEL, CADD, SIFT, PolyPhen-2, GERP, gnomAD AF, BLOSUM62 |

## Design decisions

- **Out-of-fold probabilities are the headline numbers.** Every cohort variant was in the training data of the final models. The final models' probabilities are therefore optimistic, and they appear only in a collapsed "Final model" section. The out-of-fold probabilities come from models that never saw any variant of that gene, which matches how the model would behave on a new gene.
- **XGBoost is the model shown.** It is within 0.004 ROC-AUC of random forest, and its model files are small enough to version (random forest files are about 14 MB each).
- **No retraining and no FoldX or VEP at run time.** Only the 5,000 cohort variants are supported. Anything else gets a clear error message, because new variants would need VEP annotation and a FoldX run.
- **The ClinVar label is shown separately**, as the training label, so it isn't mistaken for a prediction.

## Tests

`tests/test_app.py` uses Streamlit's `AppTest` without network access (the structure download is stubbed). It renders:
- TP53 Y220C, a pathogenic variant, which is predicted pathogenic;
- TP53 P8L, a benign variant;
- an unsupported change (R175H, not in the cohort), which must show a clear error.

## Disclaimer

The app is a research and class demonstration and is not for clinical use.
