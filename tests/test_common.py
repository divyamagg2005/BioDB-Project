from biodb import config
from biodb.common import parse_protein_change, parse_transcript, usable


def test_parse_protein_change_missense():
    assert parse_protein_change("NM_000033.4(ABCD1):c.19C>T (p.Pro7Ser)") == ("P", 7, "S")


def test_parse_protein_change_rejects_synonymous_and_nonsense():
    assert parse_protein_change("NM_000033.4(ABCD1):c.21C>T (p.Pro7Pro)") is None
    assert parse_protein_change("NM_000033.4(ABCD1):c.19C>T (p.Gln7Ter)") is None


def test_parse_protein_change_rejects_missing_or_multiple():
    assert parse_protein_change("NM_000033.4(ABCD1):c.-5C>T") is None
    assert parse_protein_change("x (p.Pro7Ser) (p.Ala8Val)") is None
    assert parse_protein_change("") is None


def test_parse_transcript():
    assert parse_transcript("NM_000033.4(ABCD1):c.19C>T (p.Pro7Ser)") == "NM_000033.4"
    assert parse_transcript("NC_000023.11:g.153724856C>T") == ""


def test_usable():
    assert usable("A")
    assert not usable("-")
    assert not usable(None)
    assert not usable("  ")


def test_config_paths_are_inside_project():
    assert config.PROCESSED_DIR.is_relative_to(config.PROJECT_ROOT)
    assert (config.PROJECT_ROOT / "pyproject.toml").exists()
