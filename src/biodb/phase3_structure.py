"""Phase 3: AlphaFold structures, residue mapping and structural features.

For each Phase 1 gene the canonical full-length AlphaFold DB model (mmCIF) is
downloaded to ``data/raw/alphafold/`` and its sequence is checked against the
UniProt sequence used in Phase 1. Every cohort variant is then mapped to its
residue, the wild-type residue is confirmed, and per-residue features are
computed:

- ``plddt``: AlphaFold confidence of the residue (0-100).
- ``plddt_window``: mean pLDDT over residues i-5..i+5.
- ``rsa``: relative solvent accessibility, residue SASA (Shrake-Rupley) over
  the theoretical maximum of Tien et al. 2013.
- ``neighbor_count``: residues whose CB (CA for glycine) lies within 10 A of
  this residue's CB, a measure of burial.
- ``secondary_structure``: helix, strand or coil, from the DSSP assignment
  stored in the AlphaFold mmCIF.

The structure files are not versioned. Their SHA-256 checksums are recorded in
``structures.csv`` and they can be re-downloaded with this module.

Run with ``python -m biodb.phase3_structure``.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

from biodb import config
from biodb.common import sha256_file, write_json

# Theoretical maximum accessible surface area (A^2), Tien et al. 2013, PLoS ONE 8:e80635.
MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0, "Q": 225.0, "E": 223.0,
    "G": 104.0, "H": 224.0, "I": 197.0, "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0,
    "P": 159.0, "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}
HELIX_TYPES = {"HELX_RH_AL_P", "HELX_RH_3T_P", "HELX_RH_PI_P"}
STRAND_TYPES = {"STRN"}
NEIGHBOR_RADIUS = 10.0
PLDDT_WINDOW = 5

STRUCTURE_FIELDS = [
    "gene", "uniprot_accession", "alphafold_entry", "alphafold_version", "structure_file",
    "structure_sha256", "residues", "mean_plddt", "sequence_matches_uniprot", "status",
]
MAPPING_FIELDS = [
    "variant_id", "gene", "uniprot_accession", "protein_position", "wt_aa", "mut_aa",
    "alphafold_entry", "structure_file", "mapping_status", "failure_reason",
    "plddt", "plddt_window", "sasa", "rsa", "neighbor_count", "secondary_structure",
]


def download(url: str, path: Path, session: requests.Session) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    for attempt in range(4):
        try:
            response = session.get(url, timeout=120)
            response.raise_for_status()
            path.parent.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix(".part")
            tmp.write_bytes(response.content)
            tmp.replace(path)
            time.sleep(0.2)
            return
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def secondary_structure(cif: dict[str, Any], length: int) -> list[str]:
    """Per-residue helix/strand/coil from the mmCIF _struct_conf records (1-based list index - 1)."""
    states = ["coil"] * length
    kinds = cif.get("_struct_conf.conf_type_id", [])
    begins = cif.get("_struct_conf.beg_label_seq_id", [])
    ends = cif.get("_struct_conf.end_label_seq_id", [])
    for kind, begin, end in zip(kinds, begins, ends):
        state = "helix" if kind in HELIX_TYPES else "strand" if kind in STRAND_TYPES else None
        if state is None:
            continue
        for position in range(int(begin), int(end) + 1):
            if 1 <= position <= length:
                states[position - 1] = state
    return states


def residue_features(cif_path: Path) -> tuple[str, dict[int, dict[str, Any]]]:
    """Sequence and per-residue features of a single-chain AlphaFold model."""
    structure = MMCIFParser(QUIET=True).get_structure("model", str(cif_path))
    chain = next(structure[0].get_chains())
    residues = [r for r in chain if r.id[0] == " "]
    sequence = "".join(seq1(r.get_resname()) for r in residues)

    ShrakeRupley().compute(structure[0], level="R")
    beta_atoms = [r["CB"] if "CB" in r else r["CA"] for r in residues]
    search = NeighborSearch(beta_atoms)
    plddt = [r["CA"].get_bfactor() for r in residues]
    ss = secondary_structure(MMCIF2Dict(str(cif_path)), len(residues))

    features: dict[int, dict[str, Any]] = {}
    for index, residue in enumerate(residues):
        number = residue.id[1]
        aa = sequence[index]
        window = plddt[max(0, index - PLDDT_WINDOW): index + PLDDT_WINDOW + 1]
        neighbors = search.search(beta_atoms[index].coord, NEIGHBOR_RADIUS, level="R")
        features[number] = {
            "aa": aa,
            "plddt": round(plddt[index], 2),
            "plddt_window": round(sum(window) / len(window), 2),
            "sasa": round(residue.sasa, 2),
            "rsa": round(residue.sasa / MAX_ASA[aa], 4) if aa in MAX_ASA else None,
            "neighbor_count": len(neighbors) - 1,
            "secondary_structure": ss[index],
        }
    return sequence, features


def validate(mapping: list[dict[str, Any]], structures: list[dict[str, Any]], cohort: list[dict[str, str]], genes: int) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []

    def check(name: str, passed: bool, detail: Any = "") -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": detail})

    ok_structures = [s for s in structures if s["status"] == "ok"]
    mapped = [m for m in mapping if m["mapping_status"] == "mapped"]
    ids = [m["variant_id"] for m in mapping]
    check("structure_per_gene", len(ok_structures) == genes, len(ok_structures))
    check("structure_sequences_match_uniprot", all(s["sequence_matches_uniprot"] for s in structures))
    check("one_row_per_cohort_variant", sorted(ids) == sorted(r["variant_id"] for r in cohort), len(ids))
    check("no_duplicate_rows", len(ids) == len(set(ids)))
    check("all_variants_mapped", len(mapped) == len(mapping),
          dict(Counter(m["failure_reason"] for m in mapping if m["failure_reason"])))
    check("plddt_in_range", all(0 <= m["plddt"] <= 100 for m in mapped))
    check("rsa_in_range", all(m["rsa"] is not None and 0 <= m["rsa"] <= 1.5 for m in mapped))
    check("neighbor_count_positive", all(m["neighbor_count"] >= 1 for m in mapped))
    return checks


def run(
    cohort_path: Path = config.PHASE1_DIR / "cohort.csv",
    genes_path: Path = config.PHASE1_DIR / "selected_genes.csv",
    uniprot_cache: Path = config.PHASE1_CACHE_DIR,
    structure_dir: Path = config.PHASE3_STRUCTURE_DIR,
    output_dir: Path = config.PHASE3_DIR,
) -> dict[str, Any]:
    from biodb.phase1_cohort import ApiClient, uniprot_entry

    with cohort_path.open(encoding="utf-8", newline="") as handle:
        cohort = list(csv.DictReader(handle))
    with genes_path.open(encoding="utf-8", newline="") as handle:
        genes = list(csv.DictReader(handle))

    client = ApiClient(cache_dir=uniprot_cache, offline=True)
    session = requests.Session()
    structures: list[dict[str, Any]] = []
    features_by_gene: dict[str, dict[int, dict[str, Any]]] = {}
    for index, gene in enumerate(genes, start=1):
        name = gene["gene"]
        url = gene["alphafold_cif_url"]
        path = structure_dir / url.rsplit("/", 1)[-1]
        download(url, path, session)
        uniprot, _ = uniprot_entry(client, name)
        sequence, features = residue_features(path)
        matches = uniprot is not None and sequence == uniprot["sequence"]
        structures.append({
            "gene": name,
            "uniprot_accession": gene["uniprot_accession"],
            "alphafold_entry": gene["alphafold_entry"],
            "alphafold_version": gene["alphafold_version"],
            "structure_file": path.relative_to(config.PROJECT_ROOT).as_posix(),
            "structure_sha256": sha256_file(path),
            "residues": len(sequence),
            "mean_plddt": round(sum(f["plddt"] for f in features.values()) / len(features), 2),
            "sequence_matches_uniprot": matches,
            "status": "ok" if matches else "sequence_mismatch",
        })
        if matches:
            features_by_gene[name] = features
        print(f"  {index}/{len(genes)} {name}: {len(sequence)} residues", flush=True)

    files = {s["gene"]: s for s in structures}
    mapping: list[dict[str, Any]] = []
    for row in cohort:
        name = row["gene"]
        position = int(row["protein_position"])
        base = {
            "variant_id": row["variant_id"], "gene": name, "uniprot_accession": row["uniprot_accession"],
            "protein_position": position, "wt_aa": row["ref_aa"], "mut_aa": row["alt_aa"],
            "alphafold_entry": files[name]["alphafold_entry"], "structure_file": files[name]["structure_file"],
        }
        residue = features_by_gene.get(name, {}).get(position)
        if name not in features_by_gene:
            reason = "structure_unusable"
        elif residue is None:
            reason = "position_outside_structure"
        elif residue["aa"] != row["ref_aa"]:
            reason = f"wt_mismatch:{residue['aa']}"
        else:
            reason = ""
        if reason:
            mapping.append({**base, "mapping_status": "failed", "failure_reason": reason})
            continue
        mapping.append({
            **base, "mapping_status": "mapped", "failure_reason": "",
            **{k: residue[k] for k in ("plddt", "plddt_window", "sasa", "rsa", "neighbor_count", "secondary_structure")},
        })

    checks = validate(mapping, structures, cohort, len(genes))
    passed = all(c["passed"] for c in checks)

    output_dir.mkdir(parents=True, exist_ok=True)
    for path, rows, fields in (
        (output_dir / "structures.csv", structures, STRUCTURE_FIELDS),
        (output_dir / "variant_structure_features.csv", mapping, MAPPING_FIELDS),
    ):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows({k: ("" if v is None else v) for k, v in r.items()} for r in rows)

    mapped = [m for m in mapping if m["mapping_status"] == "mapped"]
    labels = {r["variant_id"]: r["label"] for r in cohort}

    def summary(label: str) -> dict[str, Any]:
        rows = [m for m in mapped if labels[m["variant_id"]] == label]
        median = lambda values: sorted(values)[len(values) // 2]
        return {
            "variants": len(rows),
            "median_plddt": median([m["plddt"] for m in rows]),
            "median_rsa": median([m["rsa"] for m in rows]),
            "median_neighbor_count": median([m["neighbor_count"] for m in rows]),
            "plddt_below_70": sum(m["plddt"] < 70 for m in rows),
            "secondary_structure": dict(Counter(m["secondary_structure"] for m in rows)),
        }

    report = {
        "phase": 3,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "AlphaFold Protein Structure Database (canonical full-length models)",
        "input_cohort_sha256": sha256_file(cohort_path),
        "structures": len(structures),
        "structures_ok": sum(s["status"] == "ok" for s in structures),
        "structure_storage": structure_dir.relative_to(config.PROJECT_ROOT).as_posix() + " (not versioned; checksums in structures.csv)",
        "variants": len(mapping),
        "mapped": len(mapped),
        "failed": len(mapping) - len(mapped),
        "parameters": {"neighbor_radius_angstrom": NEIGHBOR_RADIUS, "plddt_window": PLDDT_WINDOW,
                       "max_asa": "Tien et al. 2013 theoretical", "secondary_structure": "DSSP records in AlphaFold mmCIF"},
        "by_label": {"pathogenic": summary("1"), "benign": summary("0")},
    }
    write_json(output_dir / "structure_report.json", report)
    write_json(output_dir / "validation_report.json", {"passed": passed, "checks": checks})

    print(f"\nMapped {len(mapped):,}/{len(mapping):,} variants on {report['structures_ok']} structures")
    for c in checks:
        print(f"  [{'PASS' if c['passed'] else 'FAIL'}] {c['check']}  {c['detail'] if not c['passed'] else ''}")
    return {"passed": passed, "report": report}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.parse_args(argv)
    return 0 if run()["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
