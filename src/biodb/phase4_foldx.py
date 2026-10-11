"""Phase 4: FoldX ΔΔG for every cohort variant.

Workflow per protein:

1. ``prepare``: convert the AlphaFold mmCIF from Phase 3 to PDB (chain A,
   UniProt numbering, pLDDT kept in the B-factor column).
2. ``repair``: run FoldX ``RepairPDB`` once and reuse the repaired model for
   every mutation of that protein.
3. ``mutate``: run FoldX ``BuildModel`` once per variant, in its own working
   directory, so jobs spread evenly over all cores and finished jobs are
   skipped on a rerun. The mutation string is ``<wt><chain><position><mut>``,
   e.g. ``RA175H``.
4. ``collect``: read ΔΔG (the "total energy" column of the ``Average_`` file,
   mutant minus wild type, kcal/mol) and the absolute wild-type and mutant
   energies (``Raw_`` file) by column name, and write one row per cohort
   variant, including failures with their reason.

Positive ΔΔG means the mutation is predicted to destabilise the protein. It is
an estimate from a static predicted structure, not a measurement, and it is
unreliable where AlphaFold confidence is low (pLDDT < 70); such rows are
flagged.

The FoldX binary is not part of the repository. Point ``--foldx`` (or the
``FOLDX_BIN`` environment variable) at the executable; its ``molecules``
directory must sit next to it.

Run with ``python -m biodb.phase4_foldx {prepare,repair,mutate,collect,all}``.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from biodb import config
from biodb.common import sha256_file, write_json

CHAIN = "A"
NUMBER_OF_RUNS = 3
REPAIR_TIMEOUT = 8 * 3600
MUTATE_TIMEOUT = 2 * 3600
LOW_PLDDT = 70.0

RESULT_FIELDS = [
    "variant_id", "gene", "uniprot_accession", "protein_position", "wt_aa", "mut_aa", "foldx_mutation",
    "ddg", "ddg_sd", "wt_energy", "mutant_energy", "number_of_runs", "plddt", "low_plddt",
    "status", "error", "runtime_seconds",
]


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

def work_dir() -> Path:
    return config.PHASE4_WORK_DIR


def pdb_path(gene: str) -> Path:
    return work_dir() / "structures" / f"{gene}.pdb"


def repaired_path(gene: str) -> Path:
    return work_dir() / "structures" / f"{gene}_Repair.pdb"


def job_dir(variant_id: str) -> Path:
    return work_dir() / "jobs" / variant_id.replace(":", "_").replace(".", "")


def foldx_binary(path: str | None) -> Path:
    binary = Path(path or os.environ.get("FOLDX_BIN", ""))
    if not binary.is_file():
        raise FileNotFoundError("FoldX binary not found; pass --foldx or set FOLDX_BIN")
    if not (binary.parent / "molecules").is_dir():
        raise FileNotFoundError(f"FoldX 'molecules' directory missing next to {binary}")
    return binary.resolve()


def _link_molecules(binary: Path, directory: Path) -> None:
    target = directory / "molecules"
    if not target.exists():
        try:
            target.symlink_to(binary.parent / "molecules", target_is_directory=True)
        except OSError:
            shutil.copytree(binary.parent / "molecules", target)


# ---------------------------------------------------------------------------
# FoldX output parsing
# ---------------------------------------------------------------------------

def read_fxout_table(path: Path) -> list[dict[str, str]]:
    """Rows of a FoldX .fxout table, keyed by the header line that starts with 'Pdb'."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.startswith("Pdb\t"):
            header = [h.strip() for h in line.split("\t")]
            rows = []
            for row in lines[index + 1:]:
                if not row.strip():
                    continue
                values = [v.strip() for v in row.split("\t")]
                rows.append(dict(zip(header, values)))
            return rows
    raise ValueError(f"no 'Pdb' header in {path.name}")


def parse_buildmodel(directory: Path, stem: str) -> dict[str, float]:
    """ΔΔG, its SD over runs, and mean WT/mutant total energies from BuildModel output."""
    average = read_fxout_table(directory / f"Average_{stem}.fxout")
    if len(average) != 1:
        raise ValueError(f"expected one mutant row in Average_{stem}.fxout, found {len(average)}")
    ddg = float(average[0]["total energy"])
    sd = float(average[0]["SD"])

    raw = read_fxout_table(directory / f"Raw_{stem}.fxout")
    wild = [float(r["total energy"]) for r in raw if r["Pdb"].startswith("WT_")]
    mutant = [float(r["total energy"]) for r in raw if not r["Pdb"].startswith("WT_")]
    if not wild or len(wild) != len(mutant):
        raise ValueError(f"unexpected Raw_{stem}.fxout rows: {len(wild)} WT, {len(mutant)} mutant")
    wt_energy = sum(wild) / len(wild)
    mutant_energy = sum(mutant) / len(mutant)
    if abs((mutant_energy - wt_energy) - ddg) > 0.05:
        raise ValueError(f"ΔΔG {ddg} disagrees with Raw energies {mutant_energy - wt_energy:.3f}")
    return {"ddg": ddg, "ddg_sd": sd, "wt_energy": wt_energy, "mutant_energy": mutant_energy}


def foldx_mutation(wt: str, position: int, mut: str) -> str:
    return f"{wt}{CHAIN}{position}{mut}"


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------

def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def prepare(genes: set[str] | None = None) -> None:
    """Convert each Phase 3 mmCIF to a single-chain PDB file for FoldX."""
    from Bio.PDB import MMCIFParser, PDBIO

    for structure in _read_csv(config.PHASE3_DIR / "structures.csv"):
        gene = structure["gene"]
        if genes and gene not in genes:
            continue
        out = pdb_path(gene)
        if out.exists():
            continue
        out.parent.mkdir(parents=True, exist_ok=True)
        model = MMCIFParser(QUIET=True).get_structure(gene, str(config.PROJECT_ROOT / structure["structure_file"]))
        chains = list(model[0])
        if len(chains) != 1 or chains[0].id != CHAIN:
            raise ValueError(f"{gene}: expected a single chain {CHAIN}")
        writer = PDBIO()
        writer.set_structure(model)
        writer.save(str(out))
        print(f"  prepared {gene}", flush=True)


def _repair_one(binary: Path, gene: str) -> tuple[str, str, float]:
    source = pdb_path(gene)
    directory = work_dir() / "repair" / gene
    directory.mkdir(parents=True, exist_ok=True)
    _link_molecules(binary, directory)
    shutil.copy(source, directory / source.name)
    start = time.time()
    result = subprocess.run(
        [str(binary), "--command=RepairPDB", f"--pdb={source.name}"],
        cwd=directory, capture_output=True, text=True, timeout=REPAIR_TIMEOUT,
    )
    (directory / "foldx.log").write_text(result.stdout + result.stderr, encoding="utf-8")
    produced = directory / f"{source.stem}_Repair.pdb"
    if result.returncode != 0 or not produced.exists():
        return gene, f"RepairPDB failed (exit {result.returncode})", time.time() - start
    shutil.copy(produced, repaired_path(gene))
    return gene, "ok", time.time() - start


def repair(binary: Path, genes: set[str] | None, workers: int) -> None:
    """RepairPDB once per protein, largest proteins first so they do not finish last."""
    todo = []
    for structure in _read_csv(config.PHASE3_DIR / "structures.csv"):
        gene = structure["gene"]
        if (genes and gene not in genes) or repaired_path(gene).exists():
            continue
        todo.append((int(structure["residues"]), gene))
    todo.sort(reverse=True)
    timings = work_dir() / "repair_timings.jsonl"
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_repair_one, binary, gene): gene for _, gene in todo}
        for future in as_completed(futures):
            gene, status, seconds = future.result()
            with timings.open("a", encoding="utf-8") as handle:
                handle.write(json.dumps({"gene": gene, "status": status, "seconds": round(seconds, 1)}) + "\n")
            print(f"  repair {gene}: {status} in {seconds / 60:.1f} min", flush=True)


def _mutate_one(binary: Path, variant: dict[str, str]) -> dict[str, Any]:
    directory = job_dir(variant["variant_id"])
    result_file = directory / "result.json"
    if result_file.exists():
        return json.loads(result_file.read_text(encoding="utf-8"))

    gene = variant["gene"]
    mutation = foldx_mutation(variant["wt_aa"], int(variant["protein_position"]), variant["mut_aa"])
    record: dict[str, Any] = {"variant_id": variant["variant_id"], "foldx_mutation": mutation, "number_of_runs": NUMBER_OF_RUNS}
    repaired = repaired_path(gene)
    if not repaired.exists():
        return {**record, "status": "error", "error": "no repaired structure"}

    directory.mkdir(parents=True, exist_ok=True)
    _link_molecules(binary, directory)
    shutil.copy(repaired, directory / repaired.name)
    (directory / "individual_list.txt").write_text(f"{mutation};\n", encoding="utf-8")
    start = time.time()
    try:
        run = subprocess.run(
            [str(binary), "--command=BuildModel", f"--pdb={repaired.name}",
             "--mutant-file=individual_list.txt", f"--numberOfRuns={NUMBER_OF_RUNS}"],
            cwd=directory, capture_output=True, text=True, timeout=MUTATE_TIMEOUT,
        )
        (directory / "foldx.log").write_text(run.stdout + run.stderr, encoding="utf-8")
        if run.returncode != 0:
            raise RuntimeError(f"BuildModel exit {run.returncode}")
        record.update(parse_buildmodel(directory, repaired.stem))
        record.update(status="ok", error="")
    except Exception as error:  # noqa: BLE001 - every failure is recorded per variant
        record.update(status="error", error=str(error)[:300])
    record["runtime_seconds"] = round(time.time() - start, 1)
    if record["status"] == "ok":
        # Keep only the small text outputs; the model PDBs are large.
        for path in directory.glob("*.pdb"):
            path.unlink()
        result_file.write_text(json.dumps(record) + "\n", encoding="utf-8")
    return record


def mutate(binary: Path, genes: set[str] | None, workers: int, limit_per_gene: int | None) -> None:
    variants = _read_csv(config.PHASE3_DIR / "variant_structure_features.csv")
    if genes:
        variants = [v for v in variants if v["gene"] in genes]
    if limit_per_gene:
        seen: dict[str, int] = {}
        kept = []
        for v in variants:
            seen[v["gene"]] = seen.get(v["gene"], 0) + 1
            if seen[v["gene"]] <= limit_per_gene:
                kept.append(v)
        variants = kept
    print(f"{len(variants)} BuildModel jobs on {workers} workers", flush=True)
    done = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(_mutate_one, binary, v) for v in variants]
        for future in as_completed(futures):
            done += 1
            record = future.result()
            if done % 50 == 0 or record["status"] != "ok":
                print(f"  {done}/{len(variants)} {record['variant_id']}: {record['status']} {record.get('error', '')}", flush=True)


def collect(output_dir: Path = config.PHASE4_DIR) -> dict[str, Any]:
    """One row per cohort variant, with explicit status for anything not computed."""
    variants = _read_csv(config.PHASE3_DIR / "variant_structure_features.csv")
    rows = []
    for v in variants:
        result_file = job_dir(v["variant_id"]) / "result.json"
        base = {
            "variant_id": v["variant_id"], "gene": v["gene"], "uniprot_accession": v["uniprot_accession"],
            "protein_position": v["protein_position"], "wt_aa": v["wt_aa"], "mut_aa": v["mut_aa"],
            "foldx_mutation": foldx_mutation(v["wt_aa"], int(v["protein_position"]), v["mut_aa"]),
            "plddt": v["plddt"], "low_plddt": float(v["plddt"]) < LOW_PLDDT,
        }
        if result_file.exists():
            result = json.loads(result_file.read_text(encoding="utf-8"))
            rows.append({**base, **{k: result.get(k) for k in RESULT_FIELDS if k in result}})
        else:
            rows.append({**base, "status": "missing", "error": "no completed BuildModel result"})

    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "foldx_ddg.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESULT_FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: ("" if r.get(k) is None else r.get(k)) for k in RESULT_FIELDS} for r in rows)

    ok = [r for r in rows if r.get("status") == "ok"]
    summary = {
        "phase": 4,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "foldx_command": f"BuildModel --numberOfRuns={NUMBER_OF_RUNS} on RepairPDB-repaired AlphaFold models",
        "variants": len(rows),
        "ok": len(ok),
        "not_ok": {s: sum(r.get("status") == s for r in rows) for s in ("error", "missing")},
        "ddg_sha256": sha256_file(path),
    }
    write_json(output_dir / "foldx_summary.json", summary)
    print(json.dumps(summary, indent=2))
    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("stage", choices=["prepare", "repair", "mutate", "collect", "all"])
    parser.add_argument("--foldx", help="path to the FoldX executable (default: $FOLDX_BIN)")
    parser.add_argument("--genes", help="comma-separated gene subset, e.g. for a pilot")
    parser.add_argument("--limit-per-gene", type=int, help="at most this many variants per gene (pilot)")
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args(argv)
    genes = set(args.genes.split(",")) if args.genes else None

    if args.stage in ("prepare", "all"):
        prepare(genes)
    if args.stage in ("repair", "mutate", "all"):
        binary = foldx_binary(args.foldx)
        if args.stage in ("repair", "all"):
            repair(binary, genes, args.workers)
        if args.stage in ("mutate", "all"):
            mutate(binary, genes, args.workers, args.limit_per_gene)
    if args.stage in ("collect", "all"):
        collect()
    return 0


if __name__ == "__main__":
    sys.exit(main())
