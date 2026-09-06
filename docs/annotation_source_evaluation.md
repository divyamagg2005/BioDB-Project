# Phase 5A Annotation Source Evaluation

## Cohort Audit

The frozen cohort was verified without modification:

- 2,812 variants
- 1,406 pathogenic and 1,406 benign
- BRCA1, COL2A1, GRIN2B, NSD1, TP53
- GRCh37: 1,778; GRCh38: 1,034
- All rows have usable assembly, chromosome, position, reference, and alternate fields
- 2,769 SNVs, 37 indels, and 6 inversion-labelled records

## Options

### Ensembl VEP REST

Recommended for this class project. The official POST region endpoint accepts at most 200 variants. The cohort therefore requires 9 GRCh37 batches and 6 GRCh38 batches, or 15 total requests.

VEP exposes CADD, REVEL, AlphaMissense, SIFT, PolyPhen-2, conservation, and colocated population-frequency data where available. The Phase 5B parser must use a fixed transcript policy, retain raw responses, and record missing values. The six inversion-labelled records require a preflight because standard SNV/indel assumptions may not apply.

Advantages: no large download, manageable request count, official coordinate-based service, and support for both assemblies.

Risks: rate limits, transient failures, transcript-specific duplicate scores, assembly differences, and source-specific licensing for REVEL and AlphaMissense.

### dbNSFP

dbNSFP is the broadest single source and contains most requested scores and frequency/conservation fields. However, it is a genome-scale resource containing tens of millions of variants. The full archive is unnecessary for 2,812 variants. Targeted chromosome-resource extraction is possible, but requires more storage/tooling, release pinning, coordinate handling, and license review than this class project needs.

The current dbNSFP release information was checked on 2026-09-06. No dbNSFP files were downloaded.

### gnomAD

gnomAD is authoritative for allele frequency but does not replace a computational prediction source. It can be used later as an independent AF cross-check if VEP frequency coverage is insufficient.

## Recommendation

Use Ensembl VEP REST in Phase 5B. Submit assembly-specific VCF-like records in 15 batches, save raw JSON responses, and parse only:

- CADD
- REVEL
- AlphaMissense
- SIFT
- PolyPhen-2
- one documented gnomAD allele-frequency field
- one documented Ensembl conservation field

Do not include ClinVar significance, ClinVar review status, or any clinical label returned by an annotation service as model inputs. Keep them only for audit/provenance.

Phase 5A did not call VEP, download dbNSFP, download AlphaFold structures, calculate stability, or train models. The next action is the Phase 5B request-manifest and annotation implementation.
