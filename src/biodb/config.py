"""Project-wide configuration: paths, seed, and cohort parameters.

Every phase reads its locations and shared constants from here so that no
script hard-codes an absolute path.
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
FINAL_DIR = DATA_DIR / "final"
RESULTS_DIR = PROJECT_ROOT / "results"
MODELS_DIR = PROJECT_ROOT / "models"

# Global random seed used by every sampling or modelling step.
SEED = 42

# Single reference assembly for the whole project. ClinVar's variant_summary
# lists each variant once per assembly, so mixing assemblies duplicates rows.
ASSEMBLY = "GRCh38"

# ClinVar source.
CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
CLINVAR_FILE = RAW_DIR / "clinvar" / "variant_summary.txt.gz"

# Phase 1 cohort parameters.
TARGET_GENES = 100
VARIANTS_PER_CLASS = 25
# AlphaFold DB serves a single full-length model only for proteins up to this
# length; longer human proteins are split into fragments.
MAX_PROTEIN_LENGTH = 2700

PHASE1_DIR = PROCESSED_DIR / "phase1"
PHASE1_CACHE_DIR = RAW_DIR / "phase1_api_cache"

PHASE2_DIR = PROCESSED_DIR / "phase2"
PHASE2_RAW_DIR = RAW_DIR / "phase2_vep"

PHASE3_DIR = PROCESSED_DIR / "phase3"
PHASE3_STRUCTURE_DIR = RAW_DIR / "alphafold"
