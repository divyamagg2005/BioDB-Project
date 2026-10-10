"""Test Phase 1 cohort construction."""

from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from phase1_100_gene_cohort import build_100_gene_cohort


def test_phase1():
    """Run Phase 1 cohort construction."""
    print("Starting Phase 1 test...\n")
    
    result = build_100_gene_cohort(
        target_genes=100,
        variants_per_class=25,
    )
    
    print("\n" + "=" * 80)
    print("TEST SUMMARY")
    print("=" * 80)
    print(f"Cohort file: {result['cohort_path']}")
    print(f"Genes file: {result['genes_path']}")
    print(f"Total variants: {result['total_variants']:,}")
    print(f"Total genes: {result['total_genes']}")
    print(f"Validation: {'PASSED' if result['validation_passed'] else 'FAILED'}")
    
    return result


if __name__ == "__main__":
    test_phase1()
