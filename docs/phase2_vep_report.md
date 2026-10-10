# Phase 2: VEP Annotation

## Result

- All 5,000 cohort variants were annotated through the Ensembl VEP REST API (GRCh38), with 0 failures.
- All 7 validation checks pass (`data/processed/phase2/validation_report.json`).

Run with `python -m biodb.phase2_vep`.
- Requests go in 25 batches of 200, each sorted by chromosome and position as VEP requires.
- Raw responses are saved in `data/raw/phase2_vep/` (about 96 MB, not versioned; checksums in `annotation_report.json`). An interrupted run resumes from them, and `--offline` re-parses them without network access.

## Method

- **Input.** Each variant is submitted at its VCF position with the VCF alleles. The ID column carries the cohort `variant_id`, and results are joined back on that ID. A result is also checked against the chromosome, position and allele string.
- **Transcript choice.** The chosen transcript is one of the cohort gene whose protein position and amino-acid change equal the ClinVar protein change. MANE Select is preferred, then the Ensembl canonical transcript. Results: MANE Select 4,998, another matching transcript 2.
- **Score fallback.** REVEL, PolyPhen-2, SIFT and AlphaMissense are sometimes absent on the chosen transcript only because the score's source predates that transcript. In that case the median over other transcripts of the gene with the same amino-acid change is used. REVEL rarely differs between transcripts: only 4 variants showed any spread above 0.01. Rows using this fallback are listed in `scores_from_other_transcript`, which covers 1,011 variants.
- **No label leakage.** VEP also returns ClinVar data (`clin_sig`, ClinVar IDs, phenotype flags) for co-located variants. These are never read.

## Coverage and medians

| Feature | Source field | Coverage | Median, pathogenic | Median, benign |
|---|---|---:|---:|---:|
| CADD (phred) | `cadd_phred` | 100% | 27.7 | 21.6 |
| REVEL | `revel` | 99.0% | 0.887 | 0.271 |
| AlphaMissense | `alphamissense.am_pathogenicity` | 99.0% | 0.989 | 0.109 |
| SIFT | `sift_score` | 100% | 0.00 | 0.06 |
| PolyPhen-2 | `polyphen_score` | 97.0% | 0.990 | 0.038 |
| GERP (Conservation plugin) | `conservation` | 99.4% | 1.89 | 0.97 |
| gnomAD AF, max of exome and genome | `colocated_variants[].frequencies[alt]` | 51.7% | — | — |

Missing values stay empty; nothing is imputed here.

## Points for the modelling phases

- **Allele-frequency missingness depends on the label.** gnomAD has a frequency for 85% of benign but only 18% of pathogenic variants. Absence from gnomAD is informative (pathogenic variants are rare), so Phase 5 should fill missing AF with 0 and add an "absent from gnomAD" indicator. ClinVar curators also use population frequency to call variants benign (ACMG BA1/BS1), so AF partly encodes the label. Report models with and without it.
- **Predictor circularity.** CADD, REVEL, AlphaMissense, SIFT and PolyPhen-2 are strong separators, and ClinVar submitters use such predictors as evidence (PP3/BP4). They form the conventional baseline, but some circularity is expected.
