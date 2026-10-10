# Phase 5B VEP Annotation Report

## Result

The frozen cohort was annotated through the official Ensembl VEP REST API without modifying the source cohort.

- Input variants: `2,812`
- Successfully submitted: `2,806`
- Successfully annotated/matched: `2,806`
- Excluded before submission: `6` inversion-labelled records
- Requests: `15` total, consisting of `9` GRCh37 batches and `6` GRCh38 batches
- Batch size: `200`
- Raw responses: `results/raw/vep/`
- Request log: `results/metrics/vep_request_log.jsonl`

The 37 indels were submitted as VEP region records using the ClinVar `Start` coordinate and exact VCF reference/alternate alleles. The six inversion-labelled records were not submitted because a standard VCF replacement record cannot preserve inversion semantics reliably. They remain in `annotated_project_cohort.csv` with explicit `excluded` status.

## Feature Coverage

Coverage is measured against all 2,812 cohort rows. Missing values remain blank/explicitly missing; no zeros or imputation were added.

| Feature | VEP field | Coverage | Missing |
|---|---|---:|---:|
| CADD | `transcript_consequences[].cadd_phred` | 98.4708% | 1.5292% |
| REVEL | `transcript_consequences[].revel` | 76.4936% | 23.5064% |
| AlphaMissense | `transcript_consequences[].alphamissense.am_pathogenicity` | 92.6743% | 7.3257% |
| SIFT | `transcript_consequences[].sift_score` | 93.7767% | 6.2233% |
| PolyPhen-2 | `transcript_consequences[].polyphen_score` | 82.8592% | 17.1408% |
| Allele frequency | Exact alternate allele in `colocated_variants[].frequencies`, using `af`, `gnomadg`, or `gnomade` | 42.1053% | 57.8947% |
| Conservation | `transcript_consequences[].conservation` | 99.1110% | 0.8890% |

Source documentation: [Ensembl VEP region REST endpoint](https://rest.ensembl.org/documentation/info/vep_region_post).

## Matching And Leakage

Each request record contains a stable local ID. Results were joined using that ID and validated against assembly, chromosome, and one-based `Start` coordinate. Alleles were submitted from the exact ClinVar VCF fields. No gene-only join was used, and annotation IDs are unique.

`variant_annotations.csv` excludes `ClinicalSignificance`, `clinical_significance_original`, `ClinSigSimple`, `ReviewStatus`, and `pathogenicity_label`. These fields remain available only in the merged cohort for provenance and are not annotation/model features.

The original cohort SHA-256 is unchanged:

`b9f71ac6a8f73a56f45045e64de891f36ab91e1ec9dbd1c1130612332874b412`

## Limitations

- VEP annotations are transcript-specific and coverage varies by feature.
- Population frequency is available for fewer variants than the computational scores.
- Six inversion-labelled records were excluded from VEP submission.
- No dbNSFP or other fallback database was introduced.
- No AlphaFold structures, stability calculations, feature selection, or ML training were performed.
