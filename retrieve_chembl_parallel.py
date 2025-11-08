#!/usr/bin/env python3
"""
Parallel ChEMBL bioactivity data retrieval for PLINDER systems.

This script is optimized for HPC environments (like USC CARC) with:
- Multiprocessing for parallel UniProt ID processing
- Batch processing capabilities
- Resume functionality (skip already processed)
- Memory-efficient chunked processing
- Configurable worker processes

Input: annotation_table_with_uniprot.parquet (from extract_uniprot_ids.ipynb)
Output: ChEMBL activities with all annotations for each system
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
from multiprocessing import Pool, Manager, cpu_count
from functools import partial
import logging
from datetime import datetime

# ChEMBL clients
activity = new_client.activity


def setup_logging(output_dir):
    """Setup logging to both file and console."""
    log_dir = Path(output_dir) / 'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"chembl_retrieval_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)


def query_chembl_activities_complete(uniprot_id, standard_types=['IC50', 'Ki', 'Kd', 'EC50'], 
                                     max_activities=1000):
    """
    Query ChEMBL for activities with ALL annotations in ONE step.
    
    This is a simplified version for parallel processing - no verbose output.
    
    Parameters:
    -----------
    uniprot_id : str
        UniProt accession
    standard_types : list
        Activity types to retrieve
    max_activities : int
        Maximum activities to retrieve
    
    Returns:
    --------
    pd.DataFrame or None : Complete activity data with all available annotations
    """
    try:
        # Query activities - get ALL available fields
        activities = activity.filter(
            target_organism='Homo sapiens',
            standard_type__in=standard_types,
            target_components__accession=uniprot_id
        )[:max_activities]
        
        if not activities:
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
        return df
    
    except Exception as e:
        logging.error(f"Error querying {uniprot_id}: {e}")
        return None


def process_system_group(args):
    """
    Process a group of systems with the same UniProt ID.
    
    This is the worker function for parallel processing.
    
    Parameters:
    -----------
    args : tuple
        (uniprot_id, systems_df, output_dir, standard_types, max_activities, counter, total)
    
    Returns:
    --------
    dict : Results for this UniProt ID group
    """
    uniprot_id, systems_df, output_dir, standard_types, max_activities, counter, total = args
    
    # Get process-safe logger
    logger = logging.getLogger(__name__)
    
    try:
        # Check if already processed
        output_file = Path(output_dir) / f"uniprot_{uniprot_id}" / f"{uniprot_id}_chembl_activities.parquet"
        if output_file.exists():
            counter.value += 1
            logger.info(f"[{counter.value}/{total}] Skipping {uniprot_id} - already processed")
            return {
                'uniprot_id': uniprot_id,
                'success': True,
                'skipped': True,
                'num_systems': len(systems_df),
                'reason': 'already_processed'
            }
        
        # Query ChEMBL
        logger.info(f"[{counter.value}/{total}] Querying ChEMBL for {uniprot_id} ({len(systems_df)} systems)")
        
        activities_df = query_chembl_activities_complete(
            uniprot_id,
            standard_types=standard_types,
            max_activities=max_activities
        )
        
        if activities_df is None or len(activities_df) == 0:
            counter.value += 1
            logger.warning(f"[{counter.value}/{total}] No activities found for {uniprot_id}")
            return {
                'uniprot_id': uniprot_id,
                'success': False,
                'num_systems': len(systems_df),
                'reason': 'no_activities'
            }
        
        # Add system information
        # For each activity, we'll later join with systems that use this UniProt ID
        
        # Save to parquet
        output_path = Path(output_dir) / f"uniprot_{uniprot_id}"
        output_path.mkdir(parents=True, exist_ok=True)
        
        activities_df.to_parquet(output_file, index=False)
        
        # Save metadata
        metadata = {
            'uniprot_id': uniprot_id,
            'num_systems': len(systems_df),
            'system_ids': systems_df['system_id'].tolist(),
            'total_activities': len(activities_df),
            'unique_compounds': int(activities_df['molecule_chembl_id'].nunique()),
            'unique_assays': int(activities_df['assay_chembl_id'].nunique()),
            'standard_types': activities_df['standard_type'].value_counts().to_dict(),
            'total_columns': len(activities_df.columns),
            'date_retrieved': pd.Timestamp.now().isoformat()
        }
        
        metadata_file = output_path / f"{uniprot_id}_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        counter.value += 1
        logger.info(f"[{counter.value}/{total}] ✓ {uniprot_id}: {len(activities_df)} activities, {len(systems_df)} systems")
        
        return {
            'uniprot_id': uniprot_id,
            'success': True,
            'num_systems': len(systems_df),
            'num_activities': len(activities_df),
            'num_compounds': int(activities_df['molecule_chembl_id'].nunique()),
            'num_assays': int(activities_df['assay_chembl_id'].nunique())
        }
    
    except Exception as e:
        counter.value += 1
        logger.error(f"[{counter.value}/{total}] Error processing {uniprot_id}: {e}")
        return {
            'uniprot_id': uniprot_id,
            'success': False,
            'num_systems': len(systems_df),
            'reason': str(e)
        }


def merge_activities_with_systems(annotation_file, output_dir, final_output):
    """
    Merge ChEMBL activities with PLINDER system information.
    
    This creates the final output where each activity is linked to its PLINDER systems.
    """
    logger = logging.getLogger(__name__)
    logger.info("Merging ChEMBL activities with PLINDER systems...")
    
    # Load annotation table
    df = pd.read_parquet(annotation_file)
    
    # Explode uniprot_ids (each system may have multiple UniProt IDs)
    df_exploded = df.explode('uniprot_ids')
    df_exploded = df_exploded[df_exploded['uniprot_ids'].notna()]
    
    # Collect all activities
    all_data = []
    
    for uniprot_id in df_exploded['uniprot_ids'].unique():
        activity_file = Path(output_dir) / f"uniprot_{uniprot_id}" / f"{uniprot_id}_chembl_activities.parquet"
        
        if not activity_file.exists():
            continue
        
        # Load activities
        activities = pd.read_parquet(activity_file)
        
        # Get systems for this UniProt ID
        systems = df_exploded[df_exploded['uniprot_ids'] == uniprot_id]['system_id'].unique()
        
        # Add system information
        activities['mapped_uniprot_id'] = uniprot_id
        activities['plinder_systems'] = [systems.tolist()] * len(activities)
        activities['num_plinder_systems'] = len(systems)
        
        all_data.append(activities)
    
    if not all_data:
        logger.error("No activities to merge!")
        return
    
    # Combine all
    combined = pd.concat(all_data, ignore_index=True)
    
    # Save
    output_path = Path(final_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(output_path, index=False)
    
    logger.info(f"✓ Merged data saved to {output_path}")
    logger.info(f"  Total activities: {len(combined)}")
    logger.info(f"  Unique compounds: {combined['molecule_chembl_id'].nunique()}")
    logger.info(f"  Unique UniProt IDs: {combined['mapped_uniprot_id'].nunique()}")


def main():
    parser = argparse.ArgumentParser(
        description='Parallel ChEMBL data retrieval for PLINDER systems'
    )
    parser.add_argument(
        '--input',
        default='data/annotation_table_with_uniprot.parquet',
        help='Input parquet file with UniProt mappings'
    )
    parser.add_argument(
        '--output',
        default='data/chembl_parallel',
        help='Output directory for per-UniProt results'
    )
    parser.add_argument(
        '--final-output',
        default='data/chembl_activities_merged.parquet',
        help='Final merged output file'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=None,
        help='Number of parallel workers (default: CPU count)'
    )
    parser.add_argument(
        '--chunk-size',
        type=int,
        default=None,
        help='Process in chunks of N UniProt IDs (for memory management)'
    )
    parser.add_argument(
        '--max-activities',
        type=int,
        default=1000,
        help='Maximum activities per UniProt ID'
    )
    parser.add_argument(
        '--standard-types',
        nargs='+',
        default=['IC50', 'Ki', 'Kd', 'EC50'],
        help='Activity types to retrieve'
    )
    parser.add_argument(
        '--skip-merge',
        action='store_true',
        help='Skip final merge step (only retrieve per-UniProt data)'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.output)
    
    logger.info("="*80)
    logger.info("Parallel ChEMBL Data Retrieval for PLINDER Systems")
    logger.info("="*80)
    
    # Determine number of workers
    num_workers = args.workers if args.workers else max(1, cpu_count() - 1)
    logger.info(f"Using {num_workers} parallel workers")
    
    # Load annotation table
    logger.info(f"Loading annotation table from {args.input}")
    df = pd.read_parquet(args.input)
    logger.info(f"Loaded {len(df)} systems")
    
    # Explode uniprot_ids to get one row per system-UniProt combination
    df_exploded = df.explode('uniprot_ids')
    df_exploded = df_exploded[df_exploded['uniprot_ids'].notna()]
    
    # Group by UniProt ID
    grouped = df_exploded.groupby('uniprot_ids')
    
    unique_uniprot_ids = list(grouped.groups.keys())
    logger.info(f"Found {len(unique_uniprot_ids)} unique UniProt IDs")
    
    # Prepare work items
    work_items = []
    for uniprot_id, group_df in grouped:
        work_items.append((
            uniprot_id,
            group_df.reset_index(drop=True),
            args.output,
            args.standard_types,
            args.max_activities
        ))
    
    # Process in chunks if specified
    if args.chunk_size:
        chunks = [work_items[i:i+args.chunk_size] for i in range(0, len(work_items), args.chunk_size)]
        logger.info(f"Processing in {len(chunks)} chunks of {args.chunk_size}")
    else:
        chunks = [work_items]
    
    # Process all chunks
    all_results = []
    
    for chunk_idx, chunk in enumerate(chunks, 1):
        if len(chunks) > 1:
            logger.info(f"\nProcessing chunk {chunk_idx}/{len(chunks)} ({len(chunk)} UniProt IDs)")
        
        # Create shared counter for progress tracking
        manager = Manager()
        counter = manager.Value('i', 0)
        total = len(chunk)
        
        # Add counter and total to work items
        chunk_with_counter = [item + (counter, total) for item in chunk]
        
        # Process in parallel
        with Pool(processes=num_workers) as pool:
            chunk_results = pool.map(process_system_group, chunk_with_counter)
        
        all_results.extend(chunk_results)
        
        # Brief pause between chunks
        if chunk_idx < len(chunks):
            time.sleep(5)
    
    # Summary
    logger.info("\n" + "="*80)
    logger.info("RETRIEVAL SUMMARY")
    logger.info("="*80)
    
    success_count = sum(1 for r in all_results if r['success'] and not r.get('skipped', False))
    skipped_count = sum(1 for r in all_results if r.get('skipped', False))
    failed_count = sum(1 for r in all_results if not r['success'])
    
    logger.info(f"Total UniProt IDs: {len(all_results)}")
    logger.info(f"Successfully retrieved: {success_count}")
    logger.info(f"Skipped (already processed): {skipped_count}")
    logger.info(f"Failed: {failed_count}")
    
    if success_count > 0:
        total_activities = sum(r.get('num_activities', 0) for r in all_results if r['success'])
        total_systems = sum(r.get('num_systems', 0) for r in all_results if r['success'])
        logger.info(f"\nTotal activities retrieved: {total_activities}")
        logger.info(f"Total systems covered: {total_systems}")
    
    # Save summary
    summary_file = Path(args.output) / "retrieval_summary.json"
    summary = {
        'total_uniprot_ids': len(all_results),
        'successful': success_count,
        'skipped': skipped_count,
        'failed': failed_count,
        'results': all_results,
        'date': pd.Timestamp.now().isoformat(),
        'parameters': {
            'workers': num_workers,
            'chunk_size': args.chunk_size,
            'max_activities': args.max_activities,
            'standard_types': args.standard_types
        }
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"\n✓ Summary saved to {summary_file}")
    
    # Merge results
    if not args.skip_merge:
        logger.info("\n" + "="*80)
        logger.info("MERGING RESULTS")
        logger.info("="*80)
        merge_activities_with_systems(args.input, args.output, args.final_output)
    
    logger.info("\n✓ All done!")


if __name__ == '__main__':
    main()
