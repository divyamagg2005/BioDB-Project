# Phase 7A Protein Stability Tool Evaluation

## Recommendation

Use **FoldX 5.0** as the primary tool, subject to confirming eligibility for its academic license. Use **Rosetta cartesian_ddg** only as a fallback or small sensitivity-analysis subset. Do not install either tool in Phase 7A.

The project has five validated AlphaFold mmCIF structures and 2,806 non-inversion mapped variants. The proteins range from 393 to 2,696 residues, making this a manageable small batch, but the larger multidomain proteins require stricter quality control.

## FoldX Suitability

FoldX documents `RepairPDB`, `BuildModel`, and `Stability` commands and recommends repairing a structure before further calculations. Its command-line workflow is PDB-oriented, so the downloaded AlphaFold mmCIF files should first be converted to PDB. The conversion must preserve the UniProt sequence, chain identity, and residue numbering.

FoldX can estimate single-substitution energy changes. A proposed workflow is:

1. Convert each mmCIF to PDB.
2. Verify sequence and residue numbering against the validated UniProt sequence.
3. Remove or document unsupported heteroatoms, alternate locations, waters, and nonstandard residues.
4. Run `RepairPDB` once per protein.
5. Create deterministic mutation files for each cohort variant.
6. Run `BuildModel` or the appropriate stability workflow.
7. Retain raw outputs and logs, then extract one ΔΔG row per variant.

FoldX documentation: [official manual](https://foldxsuite.crg.eu/documentation).

The software is available under academic and commercial licensing. A student project should confirm academic eligibility and accept the license terms before obtaining the executable. The binary should not be committed to the repository.

## Structure Risks

The AlphaFold files are predicted structures, not experimental structures. Expected risks include missing atoms, unsupported residues, alternate locations, chain/assembly ambiguity, flexible or disordered regions, low-confidence sites, and long multidomain proteins. FoldX recognizes only supported residue types, so preprocessing must detect anything discarded or repaired.

Low-confidence AlphaFold regions should be flagged, not silently treated as experimentally reliable. A static-model energy difference is an estimate conditioned on the selected conformation and force field.

## ΔΔG Meaning

Use the convention:

`ΔΔG = ΔG_mutant - ΔG_wild_type`

Under this convention, a positive value generally indicates a less favorable, destabilizing mutation, while a negative value generally indicates a more favorable, stabilizing mutation. The exact sign and output-column interpretation must be verified against the installed FoldX 5.0 documentation before analysis.

ΔΔG is a computational estimate, not an experimental measurement. It should be used as a protein-stability feature and never as a pathogenicity label or a substitute for the ClinVar target.

## Fallback

Rosetta cartesian_ddg is the most defensible fallback for a small validation subset, but it has heavier installation, longer runtime, and more parameters. Web predictors such as DUET/DynaMut-style services may be useful for spot checks but introduce API, coverage, licensing, and reproducibility concerns and should not be the main 2,806-variant workflow without a separate evaluation.

## Phase 7B Plan

- Confirm the FoldX academic license.
- Convert and validate the five structures.
- Run one pilot mutation per protein.
- Inspect repair, mutation, chain, residue, and output quality.
- Scale deterministically to all 2,806 mapped variants.
- Store raw outputs and explicit failure/QC statuses.
- Create a derived stability feature table without changing the frozen cohort.

Phase 7A did not install FoldX, download databases, calculate stability, or train models.
