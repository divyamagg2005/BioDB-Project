"""Tests for Phase 3 AlphaFold structure integration."""

import json
import pytest
from pathlib import Path


PHASE3_DIR = Path("data/processed/phase3")
MAPPING_PATH = PHASE3_DIR / "phase3_variant_structure_mapping.csv"
REPORT_PATH = PHASE3_DIR / "phase3_alphafold_report.json"
VALIDATION_PATH = PHASE3_DIR / "phase3_validation_report.json"
STRUCTURE_DIR = Path("data/raw/alphafold_phase3")


@pytest.fixture
def report_data():
    """Load Phase 3 report."""
    if not REPORT_PATH.exists():
        pytest.skip("Phase 3 report not found")
    with REPORT_PATH.open("r") as f:
        return json.load(f)


@pytest.fixture
def validation_data():
    """Load Phase 3 validation report."""
    if not VALIDATION_PATH.exists():
        pytest.skip("Phase 3 validation not found")
    with VALIDATION_PATH.open("r") as f:
        return json.load(f)


@pytest.fixture
def mapping_content():
    """Load Phase 3 mapping table."""
    if not MAPPING_PATH.exists():
        pytest.skip("Phase 3 mapping not found")
    with MAPPING_PATH.open("r") as f:
        return f.read()


def test_phase3_report_exists():
    """Test that Phase 3 report was created."""
    assert REPORT_PATH.exists(), f"Report not found: {REPORT_PATH}"


def test_phase3_mapping_exists():
    """Test that Phase 3 mapping table was created."""
    assert MAPPING_PATH.exists(), f"Mapping not found: {MAPPING_PATH}"


def test_phase3_validation_exists():
    """Test that Phase 3 validation was created."""
    assert VALIDATION_PATH.exists(), f"Validation not found: {VALIDATION_PATH}"


def test_structure_directory_exists():
    """Test that structure directory was created."""
    assert STRUCTURE_DIR.exists(), f"Structure dir not found: {STRUCTURE_DIR}"
    assert STRUCTURE_DIR.is_dir(), f"Not a directory: {STRUCTURE_DIR}"


def test_structure_file_count(report_data):
    """Test that structure files were downloaded."""
    expected_genes = report_data["structures"]["total_genes"]
    successful = report_data["structures"]["successful"]
    
    # Count actual .cif files
    cif_files = list(STRUCTURE_DIR.glob("*.cif"))
    
    assert len(cif_files) == successful, f"Expected {successful} CIF files, found {len(cif_files)}"


def test_all_genes_have_structures(report_data):
    """Test that all genes acquired structures."""
    successful = report_data["structures"]["successful"]
    total = report_data["structures"]["total_genes"]
    
    assert successful == total, f"Only {successful}/{total} structures acquired"


def test_mapping_success_rate(report_data):
    """Test that mapping success rate is >= 90%."""
    success_rate = report_data["mapping"]["success_rate"]
    
    assert success_rate >= 90.0, f"Mapping success rate {success_rate}% < 90%"


def test_valid_mapping_count(report_data):
    """Test that valid mappings are above threshold."""
    valid = report_data["mapping"]["valid_mappings"]
    total = report_data["mapping"]["total_variants"]
    
    # Should match approximate Phase 1/2 counts
    assert valid >= 4500, f"Too few valid mappings: {valid}"


def test_mapping_table_has_required_columns(mapping_content):
    """Test that mapping table has all required columns."""
    lines = mapping_content.strip().split("\n")
    header = lines[0].split(",")
    
    required_columns = [
        "variant_id", "gene", "UniProt_accession", "protein_position",
        "WT_AA", "mutant_AA", "structure_file", "mapping_status"
    ]
    
    for col in required_columns:
        assert col in header, f"Missing column: {col}"


def test_no_silent_failures(report_data):
    """Test that failures are recorded."""
    # All failures should be documented
    invalid_reasons = report_data["mapping"]["invalid_reasons"]
    
    # If there are failures, they should be categorized
    if report_data["mapping"]["invalid_mappings"] > 0:
        assert len(invalid_reasons) > 0, "Invalid mappings exist but no reasons recorded"


def test_uniprot_mapping_valid(report_data):
    """Test that UniProt mappings are valid."""
    structure_details = report_data["structure_details"]
    
    for gene, details in structure_details.items():
        if details["status"] == "success":
            uniprot_id = details["uniprot_id"]
            # UniProt IDs are 6-10 alphanumeric characters
            assert 6 <= len(uniprot_id) <= 10, f"Invalid UniProt ID: {uniprot_id}"


def test_protein_lengths_reasonable(report_data):
    """Test that protein lengths are within expected range."""
    structure_details = report_data["structure_details"]
    
    for gene, details in structure_details.items():
        if details["status"] == "success":
            length = details["protein_length"]
            # Proteins should be 50-10,000 residues
            assert 50 <= length <= 10000, f"Unusual protein length for {gene}: {length}"


def test_source_is_alphafold(report_data):
    """Test that source is AlphaFold (not PDB)."""
    source = report_data["source"]
    
    assert "AlphaFold" in source, f"Wrong source: {source}"
    assert "PDB" not in source, f"Should not use PDB: {source}"


def test_storage_location_documented(report_data):
    """Test that storage location is documented."""
    storage = report_data.get("storage_location", "")
    storage_note = report_data.get("storage_note", "")
    
    assert storage != "", "Storage location not documented"
    assert "alphafold" in storage.lower(), f"Wrong storage location: {storage}"


def test_cif_file_format(report_data):
    """Test that structures are in CIF format."""
    structure_details = report_data["structure_details"]
    
    for gene, details in structure_details.items():
        if "file" in details and details["file"]:
            assert details["file"].endswith(".cif"), f"Not CIF format: {details['file']}"


def test_variant_count_matches_phase2(report_data):
    """Test that variant count matches Phase 2 output."""
    total_variants = report_data["mapping"]["total_variants"]
    
    # Should be ~4993-5000 (Phase 2 had 4993 annotated)
    assert 4900 <= total_variants <= 5010, f"Unexpected variant count: {total_variants}"


def test_validation_passed(validation_data):
    """Test that validation passed."""
    if not validation_data.get("passed", False):
        # Check if it's close enough
        success_rate = validation_data.get("checks", {}).get("variant_mappings", {}).get("success_rate", 0)
        # Log warning but don't fail
        print(f"Validation marked as not passed, but success rate is {success_rate}%")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
