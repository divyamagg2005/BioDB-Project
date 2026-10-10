"""Tests for Phase 2 VEP annotation."""

import json
import pytest
from pathlib import Path


PHASE2_DIR = Path("data/processed/phase2")
COHORT_PATH = PHASE2_DIR / "phase2_annotated_cohort.csv"
REPORT_PATH = PHASE2_DIR / "phase2_annotation_report.json"


@pytest.fixture
def report_data():
    """Load Phase 2 annotation report."""
    if not REPORT_PATH.exists():
        pytest.skip("Phase 2 report not found")
    with REPORT_PATH.open("r") as f:
        return json.load(f)


@pytest.fixture
def cohort_content():
    """Load Phase 2 annotated cohort."""
    if not COHORT_PATH.exists():
        pytest.skip("Phase 2 cohort not found")
    with COHORT_PATH.open("r") as f:
        return f.read()


def test_phase2_report_exists():
    """Test that Phase 2 report was created."""
    assert REPORT_PATH.exists(), f"Report not found: {REPORT_PATH}"


def test_phase2_cohort_exists():
    """Test that Phase 2 annotated cohort was created."""
    assert COHORT_PATH.exists(), f"Cohort not found: {COHORT_PATH}"


def test_input_variants_count(report_data):
    """Test that input count matches Phase 1 (5,000 variants)."""
    input_count = report_data["input"]["variant_count"]
    assert input_count == 5000, f"Expected 5,000 variants, got {input_count}"


def test_annotation_success_rate(report_data):
    """Test that annotation success rate is >= 95%."""
    success_rate = report_data["annotation_results"]["success_rate"]
    assert success_rate >= 95.0, f"Success rate {success_rate}% < 95%"


def test_output_variant_count(report_data):
    """Test that output count is within acceptable range."""
    annotated = report_data["annotation_results"]["successfully_annotated"]
    assert annotated >= 4750, f"Too few annotations: {annotated} (expected >= 4750)"


def test_no_duplicate_annotations(report_data):
    """Test that there are no duplicate annotations in output."""
    # The validation check should confirm this
    validation = report_data["validation"]
    assert validation["checks"].get("no_duplicates", {}).get("passed", False)


def test_feature_coverage_present(report_data):
    """Test that all required features have coverage data."""
    required_features = ["CADD", "REVEL", "AlphaMissense", "SIFT", "PolyPhen2", "allele_frequency", "conservation"]
    feature_coverage = report_data["feature_coverage"]
    
    for feature in required_features:
        assert feature in feature_coverage, f"Missing feature: {feature}"


def test_core_feature_coverage_threshold(report_data):
    """Test that core features (CADD, SIFT) have >= 50% coverage."""
    feature_coverage = report_data["feature_coverage"]
    
    # These are expected to have good coverage
    for feature in ["CADD", "SIFT"]:
        coverage = feature_coverage[feature]["coverage_percentage"]
        assert coverage >= 50.0, f"{feature} coverage {coverage}% < 50%"


def test_cohort_unchanged(report_data):
    """Test that Phase 1 cohort was not modified."""
    unchanged = report_data["integrity_check"]["cohort_unchanged"]
    assert unchanged, "Phase 1 cohort was modified during annotation"


def test_no_fabricated_values(report_data):
    """Verify no fabricated values flag is set."""
    no_fabrication = report_data["integrity_check"]["no_fabricated_values"]
    assert no_fabrication, "Report indicates fabricated values"


def test_missing_values_preserved(report_data):
    """Verify missing values are preserved as null."""
    preserved = report_data["integrity_check"]["missing_values_preserved"]
    assert preserved, "Report indicates missing values were not preserved"


def test_output_cohort_has_annotation_columns(cohort_content):
    """Test that output cohort includes all annotation columns."""
    lines = cohort_content.strip().split("\n")
    header = lines[0].split(",")
    
    required_columns = ["CADD", "REVEL", "AlphaMissense", "SIFT", "PolyPhen2", "annotation_status"]
    
    for col in required_columns:
        assert col in header, f"Missing column: {col}"


def test_output_cohort_has_stable_ids(cohort_content):
    """Test that output cohort includes stable variant IDs."""
    lines = cohort_content.strip().split("\n")
    header = lines[0].split(",")
    
    assert "stable_variant_id" in header, "Missing stable_variant_id column"


def test_variant_count_approximately_equal(report_data):
    """Test input and output counts are approximately equal."""
    input_count = report_data["input"]["variant_count"]
    output_count = report_data["annotation_results"]["successfully_annotated"]
    
    # Allow up to 5% loss
    assert output_count >= input_count * 0.95, f"Too many variants lost: {input_count} -> {output_count}"


def test_batch_size_configuration(report_data):
    """Test that batching was configured correctly."""
    batch_size = report_data["requests"]["batch_size"]
    assert batch_size == 200, f"Unexpected batch size: {batch_size}"


def test_retry_configuration(report_data):
    """Test that retry logic was configured."""
    max_retries = report_data["requests"]["max_retries"]
    assert max_retries >= 3, f"Max retries too low: {max_retries}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
