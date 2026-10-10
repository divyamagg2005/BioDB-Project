# Legacy code

Archived for reference only. Nothing here is maintained or imported by `src/biodb`.

## `five_gene/`

The original 5-gene pilot (BRCA1, COL2A1, GRIN2B, NSD1, TP53; 2,812 variants), built in September 2026: ClinVar filtering, cohort selection, VEP annotation, and AlphaFold acquisition, with the reports and raw VEP responses it produced. Scripts expect to be run from the repository root with the old `data/processed/` file names.

## `hermes/`

The first attempt at the 100-gene expansion (phases 1–4, committed 2026-10-07/08). It was replaced in October 2026 because of these defects:

- Phase 1 kept both the GRCh37 and GRCh38 rows of each ClinVar variant, so the same variant was counted twice and could be sampled twice. It also applied no review-status filter and sampled the lowest-coordinate variants instead of a random subset.
- Phase 2 submitted indels at `Start` together with VCF alleles, which belong at `PositionVCF`, and picked the transcript arbitrarily.
- Phase 3 took the first AlphaFold prediction returned for each accession and did not ensure that it covers the full protein or the ClinVar transcript's isoform, so some variants failed to map.
- Phase 4 omitted the chain from FoldX mutation strings and read ΔΔG from the wrong `.fxout` column.

The replacement lives in `src/biodb/`.
