#!/usr/bin/env python3
"""
Phase 4: FoldX ΔΔG Pilot Pipeline

Prepares and runs FoldX structure repair and mutation calculations for
the BioDB-Project protein stability analysis.

This pilot operates on a small sample (50 variants) from the full dataset.

Output Schema:
- variant_id: str
- gene: str
- UniProt: str
- protein_position: int
- WT_AA: str
- mutant_AA: str
- WT_energy: float (kcal/mol)
- mutant_energy: float (kcal/mol)
- DeltaDeltaG: float (kcal/mol, mutant - WT)
- status: str (SUCCESS, ERROR, MISSING, FAILED)
- error: str (error message if failed)

FoldX 5.1 Notes:
- No rotabase.txt dependency (removed in FoldX 5)
- Requires Java 11+ runtime
- Uses BuildModel command for stability calculations
"""

import csv
import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
STRUCTURE_DIR = PROJECT_ROOT / "data/raw/alphafold_phase3"
PILOT_DIR = PROJECT_ROOT / "data/processed/phase4/pilot"
CONFIG_FILE = PILOT_DIR / "pilot_config.json"

# FoldX configuration
FOLDX_RELATIVE_PATH = Path(os.environ.get("FOLDX_DIR", "/opt/foldx"))
FOLDX_BINARY = FOLDX_RELATIVE_PATH / "foldx_20261231"


def load_pilot_config() -> dict:
    """Load pilot configuration."""
    if CONFIG_FILE.exists():
        return json.loads(CONFIG_FILE.read_text())
    
    # Default configuration
    return {
        "phase": 4,
        "pilot": {
            "enabled": True,
            "sample_size": 50,
            "genes": ["ABCD1", "ATP7A", "ATRX", "BMPR2", "BRCA1"],
            "variants_per_gene": 10
        },
        "foldx": {
            "version": "5.1",
            "binary": "foldx_20261231",
            "requires_java": True
        }
    }


def check_foldx_available() -> tuple[bool, str]:
    """Check if FoldX can run on this system."""
    # Check architecture
    machine = platform.machine()
    if machine != "x86_64":
        return False, f"FoldX requires x86_64 architecture, this system is {machine}"
    
    # Check Java
    result = subprocess.run(["java", "-version"], capture_output=True, text=True)
    if result.returncode != 0:
        return False, "Java runtime not available"
    
    # Check FoldX binary
    if not FOLDX_BINARY.exists():
        return False, f"FoldX binary not found at {FOLDX_BINARY}"
    
    # Check if executable
    if not os.access(FOLDX_BINARY, os.X_OK):
        return False, f"FoldX binary not executable: {FOLDX_BINARY}"
    
    return True, f"FoldX ready at {FOLDX_BINARY}"


def convert_cif_to_pdb(cif_file: Path, pdb_file: Path) -> tuple[bool, str]:
    """Convert AlphaFold mmCIF to PDB format for FoldX.
    
    Simple conversion that extracts ATOM records from CIF format.
    """
    atom_lines = []
    
    try:
        # CIF files may be gzipped in AlphaFold downloads
        open_func = gzip.open if cif_file.suffix == '.gz' else open
        text_mode = 'rt' if cif_file.suffix == '.gz' else 'r'
        
        with open_func(cif_file, text_mode) as f:
            in_atom_site = False
            header_fields = []
            
            for line in f:
                # Find atom site block
                if line.startswith("_atom_site."):
                    in_atom_site = True
                    # Collect field names
                    field = line.strip().split('.')[1]
                    header_fields.append(field)
                    continue
                
                if in_atom_site:
                    # End of atom site block
                    if line.startswith("_") and not line.startswith("_atom_site"):
                        in_atom_site = False
                        continue
                    
                    if line.startswith("#") or not line.strip():
                        in_atom_site = False
                        continue
                    
                    # Parse atom record
                    parts = line.split()
                    
                    if len(parts) < 20:
                        continue
                    
                    # Extract key fields (positions in CIF)
                    try:
                        # Find position of key fields
                        group_pdb_idx = header_fields.index("group_PDB") if "group_PDB" in header_fields else 0
                        atom_id_idx = header_fields.index("id") if "id" in header_fields else 1
                        label_atom_id_idx = header_fields.index("label_atom_id") if "label_atom_id" in header_fields else 3
                        label_comp_id_idx = header_fields.index("label_comp_id") if "label_comp_id" in header_fields else 5
                        auth_asym_id_idx = header_fields.index("auth_asym_id") if "auth_asym_id" in header_fields else 6
                        auth_seq_id_idx = header_fields.index("auth_seq_id") if "auth_seq_id" in header_fields else 8
                        
                        Cartn_x_idx = header_fields.index("Cartn_x") if "Cartn_x" in header_fields else 10
                        Cartn_y_idx = header_fields.index("Cartn_y") if "Cartn_y" in header_fields else 11
                        Cartn_z_idx = header_fields.index("Cartn_z") if "Cartn_z" in header_fields else 12
                        
                        # Get values
                        record_type = parts[group_pdb_idx] if group_pdb_idx < len(parts) else "ATOM"
                        atom_num = int(parts[atom_id_idx]) if atom_id_idx < len(parts) else 0
                        atom_name = parts[label_atom_id_idx] if label_atom_id_idx < len(parts) else "CA"
                        res_name = parts[label_comp_id_idx] if label_comp_id_idx < len(parts) else "ALA"
                        chain_id = parts[auth_asym_id_idx] if auth_asym_id_idx < len(parts) else "A"
                        res_num = parts[auth_seq_id_idx] if auth_seq_id_idx < len(parts) else "1"
                        
                        x = float(parts[Cartn_x_idx]) if Cartn_x_idx < len(parts) else 0.0
                        y = float(parts[Cartn_y_idx]) if Cartn_y_idx < len(parts) else 0.0
                        z = float(parts[Cartn_z_idx]) if Cartn_z_idx < len(parts) else 0.0
                        
                        # Format as PDB ATOM record
                        # ATOM format: https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
                        if record_type == "ATOM":
                            pdb_line = (
                                f"ATOM  {atom_num:>5} {atom_name:<4}{res_name:>3} "
                                f"{chain_id}{int(res_num):>4}    "
                                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
                                f"{1.0:>6.2f}{0.0:>6.2f}          {atom_name[0]:>2}\n"
                            )
                            atom_lines.append(pdb_line)
                    
                    except (ValueError, IndexError) as e:
                        # Skip malformed lines
                        continue
        
        # Write PDB file
        pdb_file.parent.mkdir(parents=True, exist_ok=True)
        with open(pdb_file, 'w') as f:
            f.writelines(atom_lines)
            f.write("TER\n")
            f.write("END\n")
        
        return True, f"Converted {len(atom_lines)} atoms"
    
    except Exception as e:
        return False, f"CIF parsing error: {e}"


def repair_structure(foldx_bin: Path, pdb_file: Path, work_dir: Path) -> tuple[bool, str, Path | None]:
    """Run FoldX RepairPDB on a structure.
    
    Returns:
        (success, message, repaired_pdb_path)
    """
    gene_name = pdb_file.stem
    repaired_file = work_dir / f"{gene_name}_Repaired.pdb"
    
    # Check if already repaired
    if repaired_file.exists():
        return True, "already_repaired", repaired_file
    
    try:
        # Run FoldX RepairPDB
        # Note: FoldX 5.1 uses different command syntax
        result = subprocess.run(
            [str(foldx_bin), "--command=RepairPDB", f"--pdb={pdb_file.name}"],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=300  # 5 minutes max per structure
        )
        
        if result.returncode == 0 and repaired_file.exists():
            return True, "repair_success", repaired_file
        else:
            return False, f"RepairPDB failed: {result.stderr[:200]}", None
    
    except subprocess.TimeoutExpired:
        return False, "RepairPDB timeout (5 min)", None
    except Exception as e:
        return False, f"RepairPDB error: {e}", None


def run_mutation(
    foldx_bin: Path,
    repaired_pdb: Path,
    mutation: str,
    work_dir: Path
) -> tuple[bool, str, float | None, float | None, float | None]:
    """Run FoldX BuildModel for a single mutation.
    
    Args:
        mutation: Format "WTposMUT" (e.g., "A100V")
    
    Returns:
        (success, message, wt_energy, mut_energy, ddg)
    """
    # Create mutation file
    mut_file = work_dir / "individual_list.txt"
    mut_file.write_text(f"{mutation};\n")
    
    try:
        # Run FoldX BuildModel
        result = subprocess.run(
            [str(foldx_bin), "--command=BuildModel", f"--pdb={repaired_pdb.name}", 
             f"--mutant-file={mut_file.name}"],
            cwd=str(work_dir),
            capture_output=True,
            text=True,
            timeout=120  # 2 minutes per mutation
        )
        
        if result.returncode != 0:
            return False, f"BuildModel failed: {result.stderr[:100]}", None, None, None
        
        # Parse FoldX output
        # FoldX writes to .fxout file
        fxout_file = work_dir / f"{repaired_pdb.stem}_1.fxout"
        
        if not fxout_file.exists():
            # Alternative output format
            fxout_files = list(work_dir.glob("*.fxout"))
            if fxout_files:
                fxout_file = fxout_files[0]
            else:
                return False, "FoldX output file not found", None, None, None
        
        # Parse energy values
        # FoldX 5.1 format: provides stability, may have different columns
        with open(fxout_file) as f:
            lines = f.readlines()
        
        # Find the mutation results
        for line in lines:
            if mutation in line:
                parts = line.split()
                # FoldX format varies; typically ΔΔG is in last column
                if len(parts) >= 2:
                    try:
                        ddg = float(parts[-1])
                        # Individual energies may be separate fields
                        wt_energy = None
                        mut_energy = None
                        if len(parts) >= 4:
                            try:
                                wt_energy = float(parts[-3])
                                mut_energy = float(parts[-2])
                            except ValueError:
                                pass
                        
                        return True, "success", wt_energy, mut_energy, ddg
                    except ValueError:
                        pass
        
        # If we didn't find the exact mutation, try parsing the last energy line
        for line in reversed(lines):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    ddg = float(parts[-1])
                    return True, "parsed_energy", None, None, ddg
                except ValueError:
                    continue
        
        return False, "Could not parse ΔΔG from output", None, None, None
    
    except subprocess.TimeoutExpired:
        return False, "BuildModel timeout (2 min)", None, None, None
    except Exception as e:
        return False, f"BuildModel error: {e}", None, None, None


def run_pilot():
    """Run the Phase 4 pilot calculation."""
    print("=" * 70)
    print("PHASE 4: FOLDX ΔΔG PILOT")
    print("=" * 70)
    print()
    
    # Load configuration
    config = load_pilot_config()
    print(f"Configuration loaded from: {CONFIG_FILE}")
    
    # Check FoldX availability
    available, message = check_foldx_available()
    print(f"\nFoldX status: {message}")
    
    if not available:
        print(f"\n⚠ Cannot run pilot locally: {message}")
        print("\nTo run on x86-64 infrastructure:")
        print("1. Deploy to AWS Fargate with x86-64 architecture")
        print("2. Or run on an x86-64 Linux machine")
        print(f"3. Set FOLDX_DIR environment variable to point to FoldX")
        print(f"   Current: {FOLDX_RELATIVE_PATH}")
        return
    
    # Load pilot sample
    pilot_file = PILOT_DIR / "pilot_sample.csv"
    if not pilot_file.exists():
        print(f"ERROR: Pilot sample not found: {pilot_file}")
        return
    
    pilot_variants = []
    with open(pilot_file) as f:
        reader = csv.DictReader(f)
        pilot_variants = list(reader)
    
    print(f"\nPilot sample: {len(pilot_variants)} variants")
    
    # Create work directory
    work_dir = PILOT_DIR / "work"
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Create results file
    results_file = PILOT_DIR / "results.csv"
    
    success_count = 0
    failed_count = 0
    
    # Track repaired structures to avoid duplication
    repaired_structures = {}
    
    print(f"\nStarting pilot run...\n")
    start_time = time.time()
    
    with open(results_file, 'w', newline='') as f:
        fieldnames = [
            'variant_id', 'gene', 'UniProt', 'protein_position',
            'WT_AA', 'mutant_AA', 'WT_energy', 'mutant_energy',
            'DeltaDeltaG', 'status', 'error'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        
        for i, variant in enumerate(pilot_variants):
            variant_id = variant['variant_id']
            gene = variant['gene']
            uniprot = variant['UniProt_accession']
            position = variant['protein_position']
            wt_aa = variant['WT_AA']
            mut_aa = variant['mutant_AA']
            
            print(f"[{i+1}/{len(pilot_variants)}] {gene} {position}{wt_aa}{mut_aa} ... ", end='', flush=True)
            
            # Get structure file
            cif_file = STRUCTURE_DIR / f"{gene}_{uniprot}.cif"
            if not cif_file.exists():
                # Try .cif.gz
                cif_file = STRUCTURE_DIR / f"{gene}_{uniprot}.cif.gz"
                if not cif_file.exists():
                    writer.writerow({
                        'variant_id': variant_id,
                        'gene': gene,
                        'UniProt': uniprot,
                        'protein_position': position,
                        'WT_AA': wt_aa,
                        'mutant_AA': mut_aa,
                        'status': 'MISSING',
                        'error': 'Structure file not found'
                    })
                    print("MISSING")
                    failed_count += 1
                    continue
            
            # Convert to PDB
            pdb_file = work_dir / f"{gene}.pdb"
            if not pdb_file.exists():
                success, msg = convert_cif_to_pdb(cif_file, pdb_file)
                if not success:
                    writer.writerow({
                        'variant_id': variant_id,
                        'gene': gene,
                        'UniProt': uniprot,
                        'protein_position': position,
                        'WT_AA': wt_aa,
                        'mutant_AA': mut_aa,
                        'status': 'CONVERSION_FAILED',
                        'error': msg
                    })
                    print(f"FAILED: {msg}")
                    failed_count += 1
                    continue
            
            # Repair structure (once per gene)
            if gene not in repaired_structures:
                repaired_pdb = work_dir / f"{gene}_Repaired.pdb"
                if not repaired_pdb.exists():
                    success, msg, repaired_file = repair_structure(FOLDX_BINARY, pdb_file, work_dir)
                    if not success:
                        writer.writerow({
                            'variant_id': variant_id,
                            'gene': gene,
                            'UniProt': uniprot,
                            'protein_position': position,
                            'WT_AA': wt_aa,
                            'mutant_AA': mut_aa,
                            'status': 'REPAIR_FAILED',
                            'error': msg
                        })
                        print(f"REPAIR FAILED: {msg}")
                        failed_count += 1
                        continue
                    else:
                        repaired_structures[gene] = repaired_file
                else:
                    repaired_structures[gene] = repaired_pdb
            
            # Run mutation
            mutation = f"{wt_aa}{position}{mut_aa}"
            success, msg, wt_energy, mut_energy, ddg = run_mutation(
                FOLDX_BINARY,
                repaired_structures[gene],
                mutation,
                work_dir
            )
            
            if success:
                writer.writerow({
                    'variant_id': variant_id,
                    'gene': gene,
                    'UniProt': uniprot,
                    'protein_position': position,
                    'WT_AA': wt_aa,
                    'mutant_AA': mut_aa,
                    'WT_energy': wt_energy if wt_energy else '',
                    'mutant_energy': mut_energy if mut_energy else '',
                    'DeltaDeltaG': ddg if ddg else '',
                    'status': 'SUCCESS',
                    'error': ''
                })
                print(f"SUCCESS (ΔΔG = {ddg:.2f} kcal/mol)")
                success_count += 1
            else:
                writer.writerow({
                    'variant_id': variant_id,
                    'gene': gene,
                    'UniProt': uniprot,
                    'protein_position': position,
                    'WT_AA': wt_aa,
                    'mutant_AA': mut_aa,
                    'status': 'FAILED',
                    'error': msg
                })
                print(f"FAILED: {msg}")
                failed_count += 1
    
    elapsed = time.time() - start_time
    
    # Summary
    print("\n" + "=" * 70)
    print("PILOT COMPLETE")
    print("=" * 70)
    print(f"  Total variants: {len(pilot_variants)}")
    print(f"  Successful: {success_count}")
    print(f"  Failed: {failed_count}")
    print(f"  Runtime: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"  Results: {results_file}")
    print("=" * 70)


if __name__ == "__main__":
    run_pilot()
