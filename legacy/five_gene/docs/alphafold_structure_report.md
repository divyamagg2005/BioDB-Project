# Phase 6 AlphaFold Structure Acquisition Report

## Source And Scope

Structures were acquired from the official [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/) API on 2026-09-07. Only the five final cohort proteins were downloaded. No unrelated structures or full-database resources were downloaded.

Format: mmCIF. API metadata and current version-6 CIF URLs were used. Sequence validation compared the downloaded CIF polymer sequence with the AlphaFold API sequence for the same reviewed human UniProt accession.

## Structures

| Gene | UniProt | AlphaFold ID | Length | File | Coverage |
|---|---|---|---:|---|---:|
| BRCA1 | P38398 | AF-P38398-F1 | 1863 | `BRCA1_P38398.cif` | 100.0% |
| COL2A1 | P02458 | AF-P02458-F1 | 1487 | `COL2A1_P02458.cif` | 100.0% |
| GRIN2B | Q13224 | AF-Q13224-F1 | 1484 | `GRIN2B_Q13224.cif` | 100.0% |
| NSD1 | Q96L73 | AF-Q96L73-F1 | 2696 | `NSD1_Q96L73.cif` | 100.0% |
| TP53 | P04637 | AF-P04637-F1 | 393 | `TP53_P04637.cif` | 100.0% |

Each accession was selected because it is the reviewed human UniProt protein used in the Phase 4A compatibility check, and the AlphaFold metadata identifies the same accession and sequence. No alternative isoform was silently substituted.

## Mapping

- Cohort rows: `2,812`
- Inversion-labelled rows excluded: `6`
- Downstream rows checked: `2,806`
- Compatible mappings: `2,806`
- Incompatible mappings: `0`
- Out-of-range positions: `0`
- Reference sequence mismatches: `0`

The six inversion-labelled records remain in the frozen cohort but are intentionally absent from `variant_structure_mapping.csv`, because their inversion semantics were excluded from downstream VEP and structure mapping.

The mapping table is [variant_structure_mapping.csv](../data/processed/variant_structure_mapping.csv). The detailed machine-readable report is [alphafold_structure_report.json](../results/metrics/alphafold_structure_report.json).

## Reproducibility

The acquisition script is [acquire_alphafold_structures.py](../src/acquire_alphafold_structures.py). It records the accession, structure ID, CIF URL, file SHA-256, API sequence length, CIF sequence length, and sequence comparison result. Variant compatibility is checked by extracting the three-letter protein change from ClinVar's protein notation and comparing the reference residue at the one-based position in the validated sequence.

No delta-G calculations, FoldX installation, stability prediction, or ML training was performed.
