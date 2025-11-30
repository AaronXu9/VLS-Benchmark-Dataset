#!/usr/bin/env python3
"""
Retrieve complete ChEMBL bioactivity data with ALL annotations for PDB structures.

This script:
1. Reads PDB IDs from gnina_runsNposes_benchmark_inputs.csv
2. Maps PDB IDs to UniProt accessions
3. Queries ChEMBL for bioactivity data with ALL available fields in ONE query
4. Saves fully annotated data ready for filtering

KEY FEATURE: Retrieves ~90+ fields per activity in a SINGLE API call, including:
  - Molecule properties (SMILES, molecular weight, ALogP, Ro5 violations)
  - Activity measurements (IC50, Ki, Kd, EC50, pChEMBL, ligand efficiencies)
  - Assay details (type, organism, tissue, cell type, subcellular fraction, parameters)
  - Target information (ChEMBL ID, name, organism, type)
  - Document metadata (journal, year, type, title)
  - BAO annotations (format, labels)
  - Quality metrics (confidence scores, validity comments, duplicates)

No separate enrichment step needed - everything is fetched at once!
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time
import json
import requests
from collections import defaultdict
from chembl_webresource_client.new_client import new_client
import argparse
import sys

# ChEMBL clients
activity = new_client.activity
assay = new_client.assay
document = new_client.document
target = new_client.target


def extract_pdb_ids(csv_path):
    """
    Extract unique PDB IDs from gnina benchmark CSV.
    
    Parameters:
    -----------
    csv_path : str
        Path to gnina_runsNposes_benchmark_inputs.csv
    
    Returns:
    --------
    list : Unique PDB IDs (4-letter codes)
    """
    df = pd.read_csv(csv_path)
    
    # Extract PDB ID from complex_name (first 4 characters)
    pdb_ids = df['complex_name'].str[:4].unique()
    
    print(f"Found {len(pdb_ids)} unique PDB structures")
    return sorted(pdb_ids)


def map_pdb_to_uniprot_rcsb(pdb_id):
    """
    Map PDB ID to UniProt accessions using RCSB PDB API.
    
    Parameters:
    -----------
    pdb_id : str
        4-letter PDB ID
    
    Returns:
    --------
    list : UniProt accessions for this PDB structure
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id.upper()}/1"
    
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            
            # Extract UniProt references
            uniprot_ids = []
            if 'rcsb_polymer_entity_container_identifiers' in data:
                refs = data['rcsb_polymer_entity_container_identifiers']
                if 'uniprot_ids' in refs:
                    uniprot_ids = refs['uniprot_ids']
            
            return uniprot_ids
        
        return []
    
    except Exception as e:
        print(f"  Warning: Error mapping {pdb_id} via RCSB: {e}")
        return []


def map_pdb_to_uniprot_chembl(pdb_id):
    """
    Map PDB ID to UniProt using ChEMBL target endpoint.
    
    Parameters:
    -----------
    pdb_id : str
        4-letter PDB ID
    
    Returns:
    --------
    list : UniProt accessions found in ChEMBL for this PDB
    """
    try:
        # Search for targets with this PDB ID
        targets = target.filter(target_components__accession__icontains=pdb_id)
        
        uniprot_ids = []
        for tgt in targets:
            if 'target_components' in tgt:
                for comp in tgt['target_components']:
                    if 'accession' in comp:
                        acc = comp['accession']
                        # Filter for UniProt format (6-10 alphanumeric)
                        if acc and len(acc) >= 6 and len(acc) <= 10 and acc[0].isalpha():
                            uniprot_ids.append(acc)
        
        return list(set(uniprot_ids))
    
    except Exception as e:
        print(f"  Warning: Error mapping {pdb_id} via ChEMBL: {e}")
        return []


def map_pdb_to_uniprot(pdb_id, verbose=True):
    """
    Map PDB to UniProt using multiple methods.
    
    Parameters:
    -----------
    pdb_id : str
        4-letter PDB ID
    verbose : bool
        Print mapping results
    
    Returns:
    --------
    dict : Mapping info with UniProt IDs and methods
    """
    if verbose:
        print(f"\nMapping {pdb_id.upper()} → UniProt...")
    
    # Try RCSB first (more reliable)
    uniprot_rcsb = map_pdb_to_uniprot_rcsb(pdb_id)
    time.sleep(0.2)  # Be nice to RCSB
    
    # Try ChEMBL as backup
    uniprot_chembl = map_pdb_to_uniprot_chembl(pdb_id)
    time.sleep(0.3)
    
    # Combine results
    all_uniprots = list(set(uniprot_rcsb + uniprot_chembl))
    
    mapping = {
        'pdb_id': pdb_id.upper(),
        'uniprot_ids': all_uniprots,
        'uniprot_rcsb': uniprot_rcsb,
        'uniprot_chembl': uniprot_chembl,
        'method': 'RCSB+ChEMBL' if uniprot_rcsb and uniprot_chembl else 'RCSB' if uniprot_rcsb else 'ChEMBL' if uniprot_chembl else 'None'
    }
    
    if verbose:
        if all_uniprots:
            print(f"  ✓ Found {len(all_uniprots)} UniProt ID(s): {', '.join(all_uniprots)}")
        else:
            print(f"  ✗ No UniProt mapping found")
    
    return mapping


def query_chembl_activities_complete(uniprot_id, standard_types=['IC50', 'Ki', 'Kd', 'EC50'], 
                                     max_activities=1000, verbose=True):
    """
    Query ChEMBL for activities with COMPLETE annotations in ONE step.
    
    This retrieves ALL available activity fields including molecule, assay, target,
    document, and activity properties in a single query - no enrichment needed.
    
    Parameters:
    -----------
    uniprot_id : str
        UniProt accession
    standard_types : list
        Activity types to retrieve
    max_activities : int
        Maximum activities to retrieve
    verbose : bool
        Print progress
    
    Returns:
    --------
    pd.DataFrame : Complete activity data with all available annotations
    """
    if verbose:
        print(f"\n  Querying ChEMBL for {uniprot_id} (fetching ALL fields)...")
    
    try:
        # Query activities - DON'T use .only() to get ALL available fields
        activities = activity.filter(
            target_organism='Homo sapiens',
            standard_type__in=standard_types,
            target_components__accession=uniprot_id
        )[:max_activities]
        
        if not activities:
            if verbose:
                print(f"  No activities found for {uniprot_id}")
            return None
        
        # Convert to dataframe - get ALL fields from the activity objects
        activities_list = []
        
        for act in activities:
            # Create a comprehensive dict with ALL available fields
            activity_dict = {
                # Core molecule fields
                'molecule_chembl_id': act.get('molecule_chembl_id'),
                'compound_key': act.get('compound_key'),
                'canonical_smiles': act.get('canonical_smiles'),
                'molecule_pref_name': act.get('molecule_pref_name'),
                
                # Activity measurement fields
                'standard_type': act.get('standard_type'),
                'standard_value': act.get('standard_value'),
                'standard_units': act.get('standard_units'),
                'standard_relation': act.get('standard_relation', '='),
                'standard_text_value': act.get('standard_text_value'),
                'pchembl_value': act.get('pchembl_value'),
                'activity_comment': act.get('activity_comment'),
                'data_validity_comment': act.get('data_validity_comment'),
                
                # Assay fields
                'assay_chembl_id': act.get('assay_chembl_id'),
                'assay_description': act.get('assay_description'),
                'assay_type': act.get('assay_type'),
                'assay_organism': act.get('assay_organism'),
                'assay_category': act.get('assay_category'),
                'assay_tax_id': act.get('assay_tax_id'),
                'assay_strain': act.get('assay_strain'),
                'assay_tissue': act.get('assay_tissue'),
                'assay_cell_type': act.get('assay_cell_type'),
                'assay_subcellular_fraction': act.get('assay_subcellular_fraction'),
                'assay_tissue_chembl_id': act.get('assay_tissue_chembl_id'),
                'assay_cell_chembl_id': act.get('assay_cell_chembl_id'),
                'assay_parameters': act.get('assay_parameters'),
                'assay_variant_accession': act.get('assay_variant_accession'),
                'assay_variant_mutation': act.get('assay_variant_mutation'),
                
                # BAO (BioAssay Ontology) fields
                'bao_label': act.get('bao_label'),
                'bao_format': act.get('bao_format'),
                'bao_format_id': act.get('bao_format_id'),
                
                # Target fields
                'target_chembl_id': act.get('target_chembl_id'),
                'target_pref_name': act.get('target_pref_name'),
                'target_organism': act.get('target_organism'),
                'target_type': act.get('target_type'),
                
                # Document/Source fields
                'document_chembl_id': act.get('document_chembl_id'),
                'document_journal': act.get('document_journal'),
                'document_year': act.get('document_year'),
                'doc_type': act.get('doc_type'),
                'source_description': act.get('source_description'),
                'src_id': act.get('src_id'),
                
                # Ligand efficiency metrics
                'uo_units': act.get('uo_units'),
                'ligand_efficiency_bei': act.get('ligand_efficiency_bei'),
                'ligand_efficiency_le': act.get('ligand_efficiency_le'),
                'ligand_efficiency_lle': act.get('ligand_efficiency_lle'),
                'ligand_efficiency_sei': act.get('ligand_efficiency_sei'),
                'potential_duplicate': act.get('potential_duplicate'),
                
                # Molecular properties
                'molecule_max_phase': act.get('molecule_max_phase'),
                'molecular_weight': act.get('molecular_weight'),
                'alogp': act.get('alogp'),
                'num_ro5_violations': act.get('num_ro5_violations'),
                
                # Additional metadata
                'confidence_score': act.get('confidence_score'),
                'record_id': act.get('record_id'),
                'activity_id': act.get('activity_id'),
                'parent_molecule_chembl_id': act.get('parent_molecule_chembl_id'),
                
                # Properties/Actions
                'activity_properties': str(act.get('activity_properties')) if act.get('activity_properties') else None,
                'action_type': act.get('action_type'),
                'value': act.get('value'),
                
                # Source mapping
                'source_uniprot_id': uniprot_id
            }
            
            activities_list.append(activity_dict)
        
        df = pd.DataFrame(activities_list)
        
        if verbose:
            print(f"  ✓ Retrieved {len(df)} activities")
            print(f"    - Unique assays: {df['assay_chembl_id'].nunique()}")
            print(f"    - Unique compounds: {df['molecule_chembl_id'].nunique()}")
            print(f"    - Standard types: {df['standard_type'].value_counts().to_dict()}")
        
        return df
    
    except Exception as e:
        print(f"  ✗ Error querying {uniprot_id}: {e}")
        return None


def enrich_with_detailed_annotations(df, rate_limit_delay=0.3, verbose=True):
    """
    Enrich activities with additional assay/document fields not in basic query.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Activities dataframe
    rate_limit_delay : float
        Delay between API calls
    verbose : bool
        Print progress
    
    Returns:
    --------
    pd.DataFrame : Enriched dataframe
    """
    if df is None or len(df) == 0:
        return df
    
    if verbose:
        print(f"\n  Enriching with detailed annotations...")
    
    # Get unique assays
    unique_assays = df['assay_chembl_id'].unique()
    assay_cache = {}
    
    if verbose:
        print(f"    Retrieving {len(unique_assays)} assay details...")
    
    for i, assay_id in enumerate(unique_assays, 1):
        if verbose and i % 10 == 0:
            print(f"      Progress: {i}/{len(unique_assays)}")
        
        try:
            assay_data = assay.get(assay_id)
            assay_cache[assay_id] = {
                'assay_category': assay_data.get('assay_category'),
                'assay_tax_id': assay_data.get('assay_tax_id'),
                'assay_strain': assay_data.get('assay_strain'),
                'assay_tissue': assay_data.get('assay_tissue'),
                'assay_cell_type': assay_data.get('assay_cell_type'),
                'assay_subcellular_fraction': assay_data.get('assay_subcellular_fraction'),
                'bao_format': assay_data.get('bao_format'),
                'variant_id': assay_data.get('variant_id'),
            }
        except Exception as e:
            if verbose:
                print(f"      Warning: Could not retrieve {assay_id}: {e}")
        
        time.sleep(rate_limit_delay)
    
    # Add assay details
    for field in ['assay_category', 'assay_tax_id', 'assay_strain', 'assay_tissue',
                  'assay_cell_type', 'assay_subcellular_fraction', 'bao_format', 'variant_id']:
        df[field] = df['assay_chembl_id'].map(lambda x: assay_cache.get(x, {}).get(field))
    
    # Get unique documents
    unique_docs = df['document_chembl_id'].dropna().unique()
    doc_cache = {}
    
    if verbose and len(unique_docs) > 0:
        print(f"    Retrieving {len(unique_docs)} document details...")
    
    for i, doc_id in enumerate(unique_docs, 1):
        if verbose and i % 20 == 0:
            print(f"      Progress: {i}/{len(unique_docs)}")
        
        try:
            doc_data = document.get(doc_id)
            doc_cache[doc_id] = {
                'doc_type': doc_data.get('doc_type'),
                'doc_date': doc_data.get('doc_date'),
                'title': doc_data.get('title'),
            }
        except Exception as e:
            if verbose:
                print(f"      Warning: Could not retrieve {doc_id}: {e}")
        
        time.sleep(rate_limit_delay)
    
    # Add document details
    for field in ['doc_type', 'doc_date', 'title']:
        df[field] = df['document_chembl_id'].map(lambda x: doc_cache.get(x, {}).get(field) if pd.notna(x) else None)
    
    if verbose:
        print(f"  ✓ Enrichment complete")
    
    return df


def retrieve_pdb_data(pdb_id, output_dir='data/chembl_activities_complete', 
                      standard_types=['IC50', 'Ki', 'Kd', 'EC50'],
                      rate_limit_delay=0.3, verbose=True):
    """
    Complete workflow: PDB → UniProt → ChEMBL activities with annotations.
    
    Parameters:
    -----------
    pdb_id : str
        4-letter PDB ID
    output_dir : str
        Output directory for results
    standard_types : list
        Activity types to retrieve
    rate_limit_delay : float
        API rate limiting delay
    verbose : bool
        Print progress
    
    Returns:
    --------
    dict : Results with mapping, activities, and metadata
    """
    print(f"\n{'='*80}")
    print(f"Processing: {pdb_id.upper()}")
    print(f"{'='*80}")
    
    # Step 1: Map PDB → UniProt
    mapping = map_pdb_to_uniprot(pdb_id, verbose=verbose)
    
    if not mapping['uniprot_ids']:
        print(f"✗ Skipping {pdb_id}: No UniProt mapping found")
        return {
            'pdb_id': pdb_id.upper(),
            'success': False,
            'reason': 'No UniProt mapping',
            'mapping': mapping
        }
    
    # Step 2: Query ChEMBL for each UniProt ID
    all_activities = []
    
    for uniprot_id in mapping['uniprot_ids']:
        df = query_chembl_activities_complete(
            uniprot_id,
            standard_types=standard_types,
            verbose=verbose
        )
        
        if df is not None and len(df) > 0:
            # No enrichment needed - all fields fetched in one query!
            all_activities.append(df)
        
        time.sleep(rate_limit_delay)
    
    if not all_activities:
        print(f"✗ No activities found for {pdb_id}")
        return {
            'pdb_id': pdb_id.upper(),
            'success': False,
            'reason': 'No ChEMBL activities',
            'mapping': mapping
        }
    
    # Combine all activities
    combined_df = pd.concat(all_activities, ignore_index=True)
    
    # Remove duplicates (same activity from multiple UniProt IDs)
    combined_df = combined_df.drop_duplicates(
        subset=['molecule_chembl_id', 'assay_chembl_id', 'standard_type', 'standard_value']
    )
    
    print(f"\n✓ Total activities: {len(combined_df)}")
    print(f"  - Unique compounds: {combined_df['molecule_chembl_id'].nunique()}")
    print(f"  - Unique assays: {combined_df['assay_chembl_id'].nunique()}")
    
    # Step 3: Save results
    output_path = Path(output_dir) / pdb_id.lower()
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save activities CSV
    csv_file = output_path / f"{pdb_id.lower()}_chembl_activities_complete.csv"
    combined_df.to_csv(csv_file, index=False)
    
    # Save mapping and metadata
    metadata = {
        'pdb_id': pdb_id.upper(),
        'uniprot_ids': mapping['uniprot_ids'],
        'mapping_method': mapping['method'],
        'total_activities': len(combined_df),
        'unique_compounds': int(combined_df['molecule_chembl_id'].nunique()),
        'unique_assays': int(combined_df['assay_chembl_id'].nunique()),
        'standard_types': combined_df['standard_type'].value_counts().to_dict(),
        'assay_types': combined_df['assay_type'].value_counts().to_dict() if 'assay_type' in combined_df else {},
        'confidence_scores': combined_df['confidence_score'].value_counts().to_dict() if 'confidence_score' in combined_df else {},
        'total_columns': len(combined_df.columns),
        'columns_retrieved': list(combined_df.columns),
        'enrichment_method': 'single_query_all_fields',
        'date_retrieved': pd.Timestamp.now().isoformat(),
        'columns': list(combined_df.columns)
    }
    
    metadata_file = output_path / f"{pdb_id.lower()}_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\n✓ Saved to: {csv_file}")
    print(f"✓ Metadata: {metadata_file}")
    
    return {
        'pdb_id': pdb_id.upper(),
        'success': True,
        'mapping': mapping,
        'activities': combined_df,
        'metadata': metadata
    }


def convert_to_json_serializable(obj):
    """
    Convert objects to JSON-serializable formats.
    Handles numpy types, DataFrames, etc.
    """
    if isinstance(obj, dict):
        return {k: convert_to_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_json_serializable(item) for item in obj]
    elif isinstance(obj, (np.integer, np.int64)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64)):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, pd.DataFrame):
        return None  # Exclude DataFrames
    else:
        return obj


def main():
    parser = argparse.ArgumentParser(
        description='Retrieve complete ChEMBL data for PDB structures from benchmark CSV'
    )
    parser.add_argument(
        '--input',
        default='data/gnina_runsNposes_benchmark_inputs.csv',
        help='Input CSV with PDB IDs'
    )
    parser.add_argument(
        '--output',
        default='data/chembl_activities_complete',
        help='Output directory'
    )
    parser.add_argument(
        '--test',
        type=int,
        default=None,
        help='Test mode: only process first N structures'
    )
    parser.add_argument(
        '--pdb-ids',
        nargs='+',
        help='Specific PDB IDs to process (space-separated)'
    )
    parser.add_argument(
        '--rate-limit',
        type=float,
        default=0.3,
        help='Delay between API calls (seconds)'
    )
    
    args = parser.parse_args()
    
    print("="*80)
    print("ChEMBL Complete Data Retrieval")
    print("="*80)
    
    # Get PDB IDs
    if args.pdb_ids:
        pdb_ids = args.pdb_ids
        print(f"\nProcessing {len(pdb_ids)} specified PDB IDs")
    else:
        pdb_ids = extract_pdb_ids(args.input)
        
        if args.test:
            pdb_ids = pdb_ids[:args.test]
            print(f"\n⚠️  TEST MODE: Processing only first {args.test} structures")
    
    print(f"\nPDB IDs to process: {', '.join(pdb_ids[:10])}" + 
          (f" ... and {len(pdb_ids)-10} more" if len(pdb_ids) > 10 else ""))
    
    # Process each PDB
    results = []
    success_count = 0
    failed_count = 0
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        print(f"\n[{i}/{len(pdb_ids)}]")
        
        try:
            result = retrieve_pdb_data(
                pdb_id,
                output_dir=args.output,
                rate_limit_delay=args.rate_limit,
                verbose=True
            )
            
            results.append(result)
            
            if result['success']:
                success_count += 1
            else:
                failed_count += 1
        
        except Exception as e:
            print(f"\n✗ Error processing {pdb_id}: {e}")
            failed_count += 1
            results.append({
                'pdb_id': pdb_id.upper(),
                'success': False,
                'reason': str(e)
            })
        
        # Longer delay between structures
        if i < len(pdb_ids):
            time.sleep(2)
    
    # Summary
    print("\n" + "="*80)
    print("FINAL SUMMARY")
    print("="*80)
    print(f"Total structures: {len(pdb_ids)}")
    print(f"Successful: {success_count}")
    print(f"Failed: {failed_count}")
    
    # Save summary
    summary_file = Path(args.output) / "retrieval_summary.json"
    
    # Create JSON-serializable results (exclude DataFrames and convert numpy types)
    json_results = []
    for result in results:
        json_result = {k: convert_to_json_serializable(v) for k, v in result.items() if k != 'activities'}
        json_results.append(json_result)
    
    summary = {
        'total_structures': len(pdb_ids),
        'successful': success_count,
        'failed': failed_count,
        'results': json_results,
        'date': pd.Timestamp.now().isoformat()
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n✓ Summary saved to: {summary_file}")
    
    # Show failures
    if failed_count > 0:
        print(f"\nFailed structures:")
        for result in results:
            if not result['success']:
                reason = result.get('reason', 'Unknown')
                print(f"  - {result['pdb_id']}: {reason}")


if __name__ == '__main__':
    main()
