# Phase 3: AlphaFold Structures and Structural Features

## Result

- 100 canonical full-length AlphaFold DB models (v6 mmCIF) were downloaded, one per cohort gene. Every model's sequence equals the UniProt sequence used in Phase 1.
- All 5,000 cohort variants map to their residue, and every wild-type residue matches.
- All 8 validation checks pass (`data/processed/phase3/validation_report.json`).

Run with `python -m biodb.phase3_structure`. Structures go to `data/raw/alphafold/`, about 126 MB, not versioned. Their SHA-256 checksums are in `structures.csv`, and rerunning the module re-downloads them.

## Features (`variant_structure_features.csv`)

| Feature | Definition |
|---|---|
| `plddt` | AlphaFold per-residue confidence (0–100), from the CA B-factor |
| `plddt_window` | Mean pLDDT over residues i−5 … i+5 |
| `sasa`, `rsa` | Residue solvent-accessible area (Biopython Shrake–Rupley), and the same divided by the theoretical maximum of Tien et al. (2013) |
| `neighbor_count` | Residues whose CB (CA for glycine) is within 10 Å of this residue's CB |
| `secondary_structure` | Helix, strand or coil, from the DSSP records in the AlphaFold mmCIF. Polyproline II, turns and bends count as coil. |

Spot check on TP53:
- R175, a known structure-destabilising hotspot, is buried (RSA 0.03, 22 neighbours).
- R248, a DNA-contact residue, is exposed (RSA 0.72).
- N- and C-terminal tail residues have pLDDT around 45.

## Pathogenic vs benign

| | Pathogenic | Benign |
|---|---:|---:|
| Median pLDDT | 89.3 | 52.4 |
| pLDDT < 70 | 446 | 1,459 |
| Median RSA | 0.14 | 0.61 |
| Median neighbour count | 19 | 8 |
| Helix / strand / coil | 1,041 / 428 / 1,031 | 563 / 217 / 1,720 |

Pathogenic variants sit mainly in confidently predicted, buried, structured regions. Benign variants sit mainly in disordered, exposed ones. Structure alone should therefore separate the classes well, and the stability-feature models must be compared against a structure baseline, not only against the VEP scores.

The 1,905 variants with pLDDT below 70 are mostly benign. FoldX ΔΔG in such regions is not meaningful, so Phase 4 should flag them rather than trust their ΔΔG.
