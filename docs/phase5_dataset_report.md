# Phase 5: Final ML Dataset

## Result

`data/processed/phase5/ml_dataset.csv` has one row per cohort variant:
- 5,000 rows across 100 genes, with 2,500 pathogenic (label 1) and 2,500 benign (label 0);
- 26 features;
- no rows lost in any merge;
- all 9 validation checks passing (`validation_report.json`).

Run with `python -m biodb.phase5_dataset`. Each feature's group and description is listed in `feature_metadata.json`. Input checksums and missingness are in `dataset_report.json`.

## Columns

- **Identifiers, not features:** `variant_id`, `gene`, `uniprot_accession`, `protein_position`, `wt_aa`, `mut_aa`. `gene` is used only to group cross-validation folds.
- **Target:** `label`.

| Group | Features |
|---|---|
| conventional (8) | `cadd_phred`, `revel`, `alphamissense`, `sift_score`, `polyphen_score`, `gerp_conservation`, `gnomad_af_max`, `gnomad_absent` |
| substitution (11) | `blosum62`, `hydropathy_change` (Kyte–Doolittle), `volume_change`, `abs_volume_change` (Zamyatnin volumes), `charge_change`, `to_/from_proline`, `to_/from_glycine`, `to_/from_cysteine` |
| structure (6) | `plddt`, `plddt_window`, `rsa`, `neighbor_count`, `ss_helix`, `ss_strand` |
| stability (1) | `ddg` (FoldX) |

## Missing values

Missing values are left empty, with no imputation in the file. Imputation happens inside the cross-validated model pipelines, fitted on training folds only.

| Feature | Missing | Pathogenic | Benign |
|---|---:|---:|---:|
| `gnomad_af_max` | 2,415 | 2,044 | 371 |
| `polyphen_score` | 151 | 69 | 82 |
| `revel` | 52 | 27 | 25 |
| `alphamissense` | 51 | 26 | 25 |
| `gerp_conservation` | 30 | 16 | 14 |

All other features are complete.

Missing gnomAD frequency is strongly tied to the label, because pathogenic alleles are rare. It is recorded explicitly as `gnomad_absent`.

## Leakage controls

- No ClinVar field (classification, review status, submitters) is a feature.
- VEP's ClinVar annotations were never extracted (Phase 2).
- Gene identity is not a feature. Models must use gene-grouped cross-validation so that no gene appears in both training and test folds.
- Remaining circularity, not removable at the data level: ClinVar curators use predictors such as REVEL/CADD (ACMG criteria PP3/BP4) and population frequency (BA1/BS1) as evidence. Phases 6–8 report models with and without these features.
