# Phase 1: 100-Gene Cohort

## Result

- 100 genes, each with 25 pathogenic and 25 benign missense variants.
- 5,000 variants in total: 2,500 pathogenic and 2,500 benign.
- All 18 validation checks pass (`data/processed/phase1/validation_report.json`).
- Rebuilding from the cached API responses (`--offline`) reproduces `cohort.csv` byte for byte. Its SHA-256 is recorded in `metadata.json`.

Run with:

```bash
python -m biodb.phase1_cohort            # first run: queries UniProt and AlphaFold DB
python -m biodb.phase1_cohort --offline  # rerun from data/raw/phase1_api_cache/
```

## Source

ClinVar `variant_summary.txt.gz`, downloaded 2026-10-10 from the NCBI FTP site. Its MD5 matched NCBI's published checksum. The file's SHA-256 is `97a20416f4d49c197a16546c8cae978794a5672306aa7130cbbd143063d80d5b`. The raw file is not versioned. A newer ClinVar release will give a different cohort, so keep this checksum with any results.

## Variant filters

Filters apply in this order; each row counts the rows that first failed at that step.

| Step | Rule | Rows removed |
|---|---|---:|
| Assembly | GRCh38 only. ClinVar lists each variant once per assembly. | 4,678,836 |
| Label | Exactly Pathogenic, Likely pathogenic, Pathogenic/Likely pathogenic (label 1) or Benign, Likely benign, Benign/Likely benign (label 0). VUS, conflicting and mixed labels are dropped. | 2,823,132 |
| Type | Single nucleotide variant | 251,876 |
| Origin | Germline | 17,164 |
| Review status | At least one star (criteria provided) | 44,303 |
| Missense | ClinVar name has exactly one `p.Xaa123Yaa` change between two different standard amino acids | 1,276,031 |
| Start loss | `p.Met1` changes are excluded because they abolish translation start rather than substitute a residue | 1,664 |
| Gene | Single gene symbol that matches the gene in the variant name | 829 |

Of 9,290,347 rows, 196,512 were kept.

Deduplication then removed 35 repeated VariationIDs and 1,324 extra records sharing a protein change within a gene. In those cases the record with the most review stars and submitters was kept. Two records whose protein change carried both labels were dropped. That left 195,151 unique missense variants.

## Gene selection

157 genes have at least 25 variants in each class. Each was mapped to exactly one reviewed human UniProt entry whose primary gene name matches. Each also needed an AlphaFold DB model of the canonical isoform covering the full protein, with a sequence identical to UniProt's.

Variants were kept only if their reference residue matches the UniProt/AlphaFold sequence at the stated position. This removed 1,624 variants, mostly ones named on a non-canonical transcript.

| Exclusion | Genes |
|---|---:|
| Protein longer than 2,700 residues (no single full-length AlphaFold model), e.g. FBN1, NF1, BRCA2, RYR1, TTN | 28 |
| Fewer than 25 per class after the residue check (isoform numbering differs), e.g. MECP2, RUNX1, TAF1 | 13 |
| Several reviewed UniProt entries for the symbol (GNAS) | 1 |

115 genes passed. They were ranked by:
1. the size of their smaller class, descending;
2. class balance, descending;
3. gene name.

The top 100 were selected. Within each gene, 25 variants per class were drawn at random with seed 42. Each gene and class has its own RNG, so the draw for one gene does not depend on which other genes were selected. The drawn variants cover the whole protein length: the median relative position is about 0.5 for both classes.

`gene_candidates.csv` records the outcome for all 157 candidate genes.

## Cohort composition

- Review stars: 1 star 3,198; 2 stars 1,494; 3 stars 308.
- Protein length: 213 to 2,696 residues (median 1,363).
- AlphaFold mean pLDDT per protein: median 68. Ten proteins have a mean below 50: BRCA1, NSD1, ANKRD11, ARID1A, CDK13, ZEB2, and the collagens COL4A1, COL4A3, COL4A4 and COL4A5. Many of their variants will lie in low-confidence regions. Per-residue pLDDT should be carried into the feature set, and ΔΔG from those regions should be treated with caution.

## Outputs (`data/processed/phase1/`)

| File | Contents |
|---|---|
| `cohort.csv` | 5,000 variants: `variant_id` (`GENE:p.R175H`), ClinVar IDs, UniProt accession, protein change, label, review status, GRCh38 VCF coordinates, ClinVar name |
| `selected_genes.csv` | 100 genes with rank, UniProt accession, length, AlphaFold entry, model version, mean pLDDT and CIF URL |
| `gene_candidates.csv` | All 157 candidate genes with mapping status and reason |
| `metadata.json` | Source checksums, parameters, and counts at every filter step |
| `validation_report.json` | The 18 acceptance checks and their results |

## Limitations

- Labels include 1-star single-submitter classifications (64% of the cohort). Requiring 2 or more stars would leave too few genes with 25 benign variants.
- Ranking by data availability favours well-studied disease genes. Results may not carry over to less-studied genes.
- ClinVar classifications often cite computational predictors (ACMG criteria PP3/BP4). Models using those predictors partly re-learn the label. This must be addressed in the modelling phases.
