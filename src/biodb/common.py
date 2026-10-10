"""Small helpers shared across phases."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

MISSING = {"", ".", "-", "NA", "N/A", "na", "n/a", "null", "None"}

AA_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}

# ClinVar names look like "NM_000033.4(ABCD1):c.19C>T (p.Pro7Ser)".
PROTEIN_CHANGE_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
TRANSCRIPT_RE = re.compile(r"^(N[MR]_\d+\.\d+)")


def usable(value: Any) -> bool:
    """True if value is present and not a missing-value token."""
    return value is not None and str(value).strip() not in MISSING


def parse_protein_change(name: str) -> tuple[str, int, str] | None:
    """Parse a ClinVar name into (ref_aa, position, alt_aa) one-letter codes.

    Returns None unless the name holds exactly one standard missense change
    with distinct reference and alternate residues.
    """
    matches = PROTEIN_CHANGE_RE.findall(name or "")
    if len(matches) != 1:
        return None
    ref3, pos, alt3 = matches[0]
    ref, alt = AA_THREE_TO_ONE.get(ref3), AA_THREE_TO_ONE.get(alt3)
    if ref is None or alt is None or ref == alt:
        return None
    return ref, int(pos), alt


def parse_transcript(name: str) -> str:
    """Return the RefSeq transcript accession from a ClinVar name, or ''."""
    match = TRANSCRIPT_RE.match(name or "")
    return match.group(1) if match else ""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=False) + "\n", encoding="utf-8")
