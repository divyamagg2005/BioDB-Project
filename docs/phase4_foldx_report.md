# Phase 4: FoldX ΔΔG

## Result

- ΔΔG was computed for all 5,000 cohort variants, with 0 errors and 0 missing.
- Output: `data/processed/phase4/foldx_ddg.csv`, SHA-256 `53bd4db8144903a0fe708757e8a51508962dc2eb43b6a5147c5b142ea2bb8af3`.

## Method

`python -m biodb.phase4_foldx all` (FoldX 5.1, binary `foldx_20261231`). It has four stages:

1. **prepare:** each Phase 3 AlphaFold mmCIF is written as a single-chain PDB (chain A, UniProt numbering, pLDDT in the B-factor column) with Biopython.
2. **repair:** `RepairPDB` runs once per protein. The repaired model is reused for every mutation of that protein.
3. **mutate:** `BuildModel --numberOfRuns=3` runs once per variant, each in its own directory. The mutation string carries the chain, e.g. `YA220C`.
4. **collect:**
   - ΔΔG is the `total energy` column of `Average_*.fxout`: the mean over 3 runs of mutant minus wild type, in kcal/mol.
   - Its standard deviation and the absolute WT and mutant energies come from `Raw_*.fxout`. Columns are read by name.
   - The parser checks that ΔΔG equals mutant minus WT energy.

Positive ΔΔG means the mutation is predicted to destabilise the protein.

**Compute:** one AWS `c7a.8xlarge` spot instance (32 vCPU, us-east-1).
- Wall time: 1 h 22 min. Total CPU time: about 42 hours.
- RepairPDB was slowest for MYH9, MTOR and ABCA12 (25–27 min each). BuildModel took a median 14 s per variant, up to 161 s.
- The instance was terminated after the results were copied.

The FoldX binary and licence are not in the repository. Raw FoldX text outputs (no PDB models) are archived locally in `data/raw/phase4_foldx/foldx_outputs.tar.gz`, which is not versioned.

## Checks

- The parser was tested on real output. TP53 Y220C gives ΔΔG +3.48 kcal/mol, close to the roughly 4 kcal/mol measured experimentally for this well-known destabilising mutant.
- Every cohort variant has exactly one row with an explicit status.
- The SD over the 3 runs is small (median 0.04 kcal/mol). A few large-clash mutations reach up to 11.8.

## Distribution

| Group | n | Median ΔΔG | Fraction > 2 kcal/mol |
|---|---:|---:|---:|
| Pathogenic, pLDDT ≥ 70 | 2,054 | 2.57 | 0.57 |
| Pathogenic, pLDDT < 70 | 446 | 0.25 | 0.22 |
| Benign, pLDDT ≥ 70 | 1,041 | 0.57 | 0.15 |
| Benign, pLDDT < 70 | 1,459 | 0.19 | 0.10 |

ΔΔG ranges from −4.2 to +45.4 kcal/mol; very large values come from buried small-to-large substitutions that FoldX scores as steric clashes. ΔΔG alone separates the classes with ROC-AUC 0.733 on the whole cohort (no model, no cross-validation).

## Caveats

- ΔΔG is computed on a static predicted structure. In low-confidence regions (pLDDT < 70, flagged in `low_plddt`) it is close to zero for both classes and carries little information.
- ΔΔG correlates with burial, so its added value has to be judged against the Phase 3 structural features as well as the VEP scores.
- The binary expires on 31 Dec 2026. Any re-run needs a newer FoldX release, which may give slightly different values.
