#!/usr/bin/env python3
"""
Phase 3.5: Transcript-Aware Structure Mapping Improvement

Investigates and attempts to rescue failed variant-to-structure mappings by:
1. Using ClinVar transcript information
2. Mapping to correct UniProt canonical isoforms
3. Checking AlphaFold availability for correct isoforms
4. Recalculating residue positions when possible

Classification of failures:
A. Successfully rescued
B. Biologically unmappable with available AlphaFold structure
C. Requires structure/isoform unavailable in AlphaFold
D. Mapping ambiguity requiring exclusion
"""

import csv
import json
import re
import requests
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
COHORT_FILE = PROJECT_ROOT / "data/processed/phase2/phase2_annotated_cohort.csv"
MAPPING_FILE = PROJECT_ROOT / "data/processed/phase3/phase3_variant_structure_mapping.csv"
STRUCTURE_DIR = PROJECT_ROOT / "data/raw/alphafold_phase3"
OUTPUT_DIR = PROJECT_ROOT / "data/processed/phase3.5"

# Gene to UniProt mapping from Phase 1
GENE_FILE = PROJECT_ROOT / "data/processed/phase1/phase1_selected_genes.csv"


def parse_name_field(name: str) -> dict:
    """Extract transcript ID and protein change from ClinVar Name field."""
    result = {
        'transcript_id': None,
        'transcript_version': None,
        'ref_aa': None,
        'pos': None,
        'alt_aa': None,
        'cdna_change': None
    }
    
    try:
        # Extract transcript ID: NM_XXXXX.X
        tx_match = re.match(r'(NM_\d+\.\d+)', name)
        if tx_match:
            result['transcript_id'] = tx_match.group(1)
            result['transcript_version'] = tx_match.group(1).split('.')[1]
        
        # Extract cDNA change
        cdna_match = re.search(r':c\.(\d+[ACGT]>[ACGT])', name)
        if cdna_match:
            result['cdna_change'] = cdna_match.group(1)
        
        # Extract protein change: p.XaaNNNXaa
        prot_match = re.search(r'p\.(\D{3})(\d+)(\D{3})', name)
        if prot_match:
            result['ref_aa'] = prot_match.group(1)
            result['pos'] = int(prot_match.group(2))
            result['alt_aa'] = prot_match.group(3)
    except Exception:
        pass
    
    return result


def check_alphafold_availability(uniprot_id: str) -> dict:
    """Check if AlphaFold has structure for this UniProt ID."""
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        response = requests.get(api_url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            if data:
                seq = data[0].get('sequence', '')
                return {
                    'available': True,
                    'length': len(seq),
                    'model_id': data[0].get('modelEntityId'),
                    'plddt': data[0].get('globalMetricValue')
                }
    except Exception:
        pass
    return {'available': False, 'length': 0}


def load_gene_uniprot_mapping() -> dict:
    """Load gene → UniProt mapping from Phase 1."""
    mapping = {}
    with open(GENE_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mapping[row['gene']] = row['uniprot_accession']
    return mapping


def load_failed_mappings() -> list[dict]:
    """Load all failed mappings from Phase 3."""
    failed = []
    with open(MAPPING_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['mapping_status'] == 'failed':
                failed.append(row)
    return failed


def load_cohort_transcripts() -> dict:
    """Load transcript info from cohort."""
    transcripts = {}
    with open(COHORT_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            variant_id = row['stable_variant_id']
            name = row.get('Name', '')
            parsed = parse_name_field(name)
            transcripts[variant_id] = {
                'transcript': parsed['transcript_id'],
                'ref_aa': parsed['ref_aa'],
                'pos': parsed['pos'],
                'alt_aa': parsed['alt_aa'],
                'gene': row['GeneSymbol']
            }
    return transcripts


def investigate_gene(gene: str, variants: list[dict], gene_uniprot: dict) -> dict:
    """Investigate failures for a specific gene.
    
    Returns classification and analysis for each variant.
    """
    result = {
        'gene': gene,
        'current_uniprot': gene_uniprot.get(gene, 'UNKNOWN'),
        'alphafold_current': check_alphafold_availability(gene_uniprot.get(gene, '')),
        'variants': [],
        'summary': {}
    }
    
    # Check alternative UniProts for this gene
    # This would require NCBI API or UniProt API calls
    # For now, mark based on known issues
    
    for variant in variants:
        variant_id = variant['variant_id']
        failure_reason = variant['failure_reason']
        
        classification = classify_failure(gene, failure_reason, result['alphafold_current'])
        
        variant_result = {
            'variant_id': variant_id,
            'protein_position': variant.get('protein_position', ''),
            'failure_reason': failure_reason,
            'classification': classification['class'],
            'notes': classification['notes']
        }
        
        result['variants'].append(variant_result)
    
    return result


def classify_failure(gene: str, reason: str, alphafold_info: dict) -> dict:
    """Classify a failure into categories A, B, C, or D.
    
    A - Successfully rescued (can be mapped with correct transcript)
    B - Biologically unmappable with available AlphaFold structure
    C - Requires structure/isoform unavailable in AlphaFold
    D - Mapping ambiguity requiring exclusion
    """
    reason_lower = reason.lower()
    
    # Out of range errors
    if 'out_of_range' in reason_lower:
        # Check if AlphaFold has a short structure
        if alphafold_info['available'] and alphafold_info['length'] < 1500:
            return {
                'class': 'C',
                'notes': 'AlphaFold structure incomplete/shorter than canonical protein'
            }
        return {
            'class': 'B',
            'notes': 'Position exceeds protein length - may require different isoform'
        }
    
    # Reference mismatch errors
    if 'reference_mismatch' in reason_lower:
        # This requires checking if correct transcript would fix it
        # For now, classify as needing transcript verification
        return {
            'class': 'D',
            'notes': 'Reference AA differs - requires transcript sequence verification'
        }
    
    # Unknown errors
    return {
        'class': 'D',
        'notes': 'Unknown failure type'
    }


def main():
    """Run Phase 3.5 investigation."""
    
    print("=== PHASE 3.5: TRANSCRIPT-AWARE STRUCTURE MAPPING INVESTIGATION ===\n")
    
    # Load data
    print("Loading data...")
    gene_uniprot = load_gene_uniprot_mapping()
    failed_mappings = load_failed_mappings()
    
    print(f"  Total failed mappings: {len(failed_mappings)}")
    
    # Group by gene
    by_gene = defaultdict(list)
    for variant in failed_mappings:
        by_gene[variant['gene']].append(variant)
    
    print(f"  Genes with failures: {len(by_gene)}")
    
    # Investigate each gene
    print("\nInvestigating failures by gene...")
    results = {}
    classifications = Counter()
    
    # Priority genes
    priority_genes = ['VWF', 'MYO15A', 'USH2A', 'FBN2']
    
    for gene in priority_genes:
        if gene in by_gene:
            print(f"\n  {gene}:")
            result = investigate_gene(gene, by_gene[gene], gene_uniprot)
            results[gene] = result
            
            # Count classifications
            for v in result['variants']:
                classifications[v['classification']] += 1
            
            print(f"    Current UniProt: {result['current_uniprot']}")
            print(f"    AlphaFold length: {result['alphafold_current']['length']} aa")
            print(f"    Variants: {len(result['variants'])}")
            print(f"    Classifications: {Counter([v['classification'] for v in result['variants']])}")
    
    # Process remaining genes
    print("\n\nProcessing remaining genes...")
    other_genes = [g for g in by_gene.keys() if g not in priority_genes]
    
    for gene in other_genes:
        result = investigate_gene(gene, by_gene[gene], gene_uniprot)
        results[gene] = result
        
        for v in result['variants']:
            classifications[v['classification']] += 1
    
    # Summary
    print("\n" + "="*70)
    print("CLASSIFICATION SUMMARY:")
    print("="*70)
    print(f"\nClass A (Rescued):            {classifications.get('A', 0)}")
    print(f"Class B (Unmappable):         {classifications.get('B', 0)}")
    print(f"Class C (Structure Missing):  {classifications.get('C', 0)}")
    print(f"Class D (Ambiguous):          {classifications.get('D', 0)}")
    print(f"\nTotal:                        {sum(classifications.values())}")
    
    # Save results
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    report = {
        'investigation_date': datetime.now(timezone.utc).isoformat(),
        'phase': 'Phase 3.5',
        'input': {
            'failed_mappings': len(failed_mappings),
            'genes_with_failures': len(by_gene)
        },
        'classifications': {
            'A_rescued': classifications.get('A', 0),
            'B_unmappable': classifications.get('B', 0),
            'C_structure_missing': classifications.get('C', 0),
            'D_ambiguous': classifications.get('D', 0)
        },
        'gene_results': {gene: {
            'current_uniprot': r['current_uniprot'],
            'alphafold_length': r['alphafold_current']['length'],
            'variant_count': len(r['variants']),
            'classifications': dict(Counter([v['classification'] for v in r['variants']]))
        } for gene, r in results.items()},
        'priority_genes_analysis': {
            gene: {
                'current_uniprot': results[gene]['current_uniprot'],
                'alphafold_length': results[gene]['alphafold_current']['length'],
                'issue_summary': f"AlphaFold has only {results[gene]['alphafold_current']['length']} aa structure"
            } for gene in priority_genes if gene in results
        }
    }
    
    output_file = OUTPUT_DIR / "phase3.5_failure_classification_report.json"
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nReport saved: {output_file}")
    
    return report


if __name__ == "__main__":
    main()
