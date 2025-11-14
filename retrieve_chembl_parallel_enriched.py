#!/usr/bin/env python3
"""
Enriched parallel ChEMBL bioactivity data retrieval for PLINDER systems.

This script uses MULTIPLE ChEMBL clients for comprehensive data:
- Activity client: bioactivity measurements (primary)
- Assay client: complete assay details (enrichment)
- Target client: target information (optional)
- Molecule client: compound properties (optional)

Key improvements over basic retrieval:
- Complete assay-level information (organism, cell type, tissue, etc.)
- Confidence scores for quality filtering
- BAO ontology annotations
- Source database information
- Links to cell lines and tissues

Input: annotation_table_with_uniprot.parquet
Output: Enriched ChEMBL activities with complete annotations
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
activity_client = new_client.activity
assay_client = new_client.assay
target_client = new_client.target
molecule_client = new_client.molecule


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


def query_chembl_activities_enriched(uniprot_id, standard_types=['IC50', 'Ki', 'Kd', 'EC50'], 
                                     max_activities=None, enrich_assays=True,
                                     enrich_targets=False, enrich_molecules=False):
    """
    Query ChEMBL for activities with enriched annotations from multiple clients.
    
    Parameters:
    -----------
    uniprot_id : str
        UniProt accession
    standard_types : list
        Activity types to retrieve
    max_activities : int or None
        Maximum activities to retrieve per UniProt ID
    enrich_assays : bool
        If True, enrich with detailed assay information (RECOMMENDED)
    enrich_targets : bool
        If True, enrich with target information
    enrich_molecules : bool
        If True, enrich with molecule properties
    
    Returns:
    --------
    pd.DataFrame or None : Complete activity data with enriched annotations
    """
    logger = logging.getLogger(__name__)
    
    try:
        # Step 1: Query activities
        activities = activity_client.filter(
            target_organism='Homo sapiens',
            standard_type__in=standard_types,
            target_components__accession=uniprot_id
        )
        
        if max_activities is not None:
            activities = activities[:max_activities]
        
        if not activities:
            return None
        
        # Convert to list of dicts
        activities_list = [dict(act) for act in activities]
        
        # Step 2: Enrich with assay details (STRONGLY RECOMMENDED)
        if enrich_assays:
            # Get unique assay IDs
            unique_assay_ids = list(set(act['assay_chembl_id'] for act in activities_list))
            logger.info(f"{uniprot_id}: Enriching {len(unique_assay_ids)} unique assays...")
            
            # Retrieve assay details
            assay_details = {}
            for assay_id in unique_assay_ids:
                try:
                    assay = assay_client.get(assay_id)
                    assay_details[assay_id] = {
                        # Assay identification
                        'assay_description_enriched': assay.get('description'),
                        'assay_type_description': assay.get('assay_type_description'),
                        
                        # Assay conditions (often missing from activity client)
                        'assay_organism_enriched': assay.get('assay_organism'),
                        'assay_strain': assay.get('assay_strain'),
                        'assay_tissue_enriched': assay.get('assay_tissue'),
                        'assay_cell_type_enriched': assay.get('assay_cell_type'),
                        'assay_subcellular_fraction_enriched': assay.get('assay_subcellular_fraction'),
                        'assay_tax_id_enriched': assay.get('assay_tax_id'),
                        'assay_test_type': assay.get('assay_test_type'),
                        
                        # Quality metrics (NOT in activity client)
                        'assay_confidence_score': assay.get('confidence_score'),
                        'assay_confidence_description': assay.get('confidence_description'),
                        'assay_relationship_type': assay.get('relationship_type'),
                        'assay_relationship_description': assay.get('relationship_description'),
                        
                        # BAO ontology
                        'assay_bao_format_enriched': assay.get('bao_format'),
                        'assay_bao_label_enriched': assay.get('bao_label'),
                        
                        # Links
                        'assay_cell_chembl_id': assay.get('cell_chembl_id'),
                        'assay_tissue_chembl_id': assay.get('tissue_chembl_id'),
                        'assay_target_chembl_id': assay.get('target_chembl_id'),
                        'assay_document_chembl_id': assay.get('document_chembl_id'),
                        
                        # Source info
                        'assay_src_id': assay.get('src_id'),
                        'assay_src_assay_id': assay.get('src_assay_id'),
                        
                        # Additional
                        'assay_category': assay.get('assay_category'),
                        'assay_parameters': str(assay.get('assay_parameters')) if assay.get('assay_parameters') else None,
                        'assay_classifications': str(assay.get('assay_classifications')) if assay.get('assay_classifications') else None,
                        'assay_variant_sequence': assay.get('variant_sequence')
                    }
                    time.sleep(0.1)  # Rate limiting
                except Exception as e:
                    logger.warning(f"Failed to enrich assay {assay_id}: {e}")
            
            # Merge assay details into activities
            for act in activities_list:
                assay_id = act['assay_chembl_id']
                if assay_id in assay_details:
                    act.update(assay_details[assay_id])
        
        # Step 3: Enrich with target details (OPTIONAL)
        if enrich_targets:
            unique_target_ids = list(set(act.get('target_chembl_id') for act in activities_list if act.get('target_chembl_id')))
            logger.info(f"{uniprot_id}: Enriching {len(unique_target_ids)} unique targets...")
            
            target_details = {}
            for target_id in unique_target_ids:
                try:
                    target = target_client.get(target_id)
                    target_details[target_id] = {
                        'target_pref_name_enriched': target.get('pref_name'),
                        'target_organism_enriched': target.get('organism'),
                        'target_type_enriched': target.get('target_type'),
                        'target_tax_id_enriched': target.get('tax_id'),
                        'target_species_group_flag': target.get('species_group_flag'),
                        'target_components_enriched': str(target.get('target_components')) if target.get('target_components') else None
                    }
                    time.sleep(0.1)
                except Exception as e:
                    logger.warning(f"Failed to enrich target {target_id}: {e}")
            
            for act in activities_list:
                target_id = act.get('target_chembl_id')
                if target_id in target_details:
                    act.update(target_details[target_id])
        
        # Step 4: Enrich with molecule details (OPTIONAL)
        if enrich_molecules:
            unique_molecule_ids = list(set(act.get('molecule_chembl_id') for act in activities_list if act.get('molecule_chembl_id')))
            logger.info(f"{uniprot_id}: Enriching {len(unique_molecule_ids)} unique molecules...")
            
            molecule_details = {}
            for mol_id in unique_molecule_ids:
                try:
                    mol = molecule_client.get(mol_id)
                    molecule_details[mol_id] = {
                        'molecule_pref_name_enriched': mol.get('pref_name'),
                        'molecule_max_phase': mol.get('max_phase'),
                        'molecule_type': mol.get('molecule_type'),
                        'molecule_first_approval': mol.get('first_approval'),
                        'molecule_therapeutic_flag': mol.get('therapeutic_flag'),
                        'molecule_natural_product': mol.get('natural_product'),
                        'molecule_properties': str(mol.get('molecule_properties')) if mol.get('molecule_properties') else None,
                        'molecule_structures': str(mol.get('molecule_structures')) if mol.get('molecule_structures') else None
                    }
                    time.sleep(0.1)
                except Exception as e:
                    logger.warning(f"Failed to enrich molecule {mol_id}: {e}")
            
            for act in activities_list:
                mol_id = act.get('molecule_chembl_id')
                if mol_id in molecule_details:
                    act.update(molecule_details[mol_id])
        
        # Convert to DataFrame
        df = pd.DataFrame(activities_list)
        
        # Add source UniProt ID
        df['source_uniprot_id'] = uniprot_id
        
        # Log statistics
        logger.info(f"{uniprot_id}: Retrieved {len(df)} activities, {df['molecule_chembl_id'].nunique()} unique compounds, {df['assay_chembl_id'].nunique()} unique assays")
        
        if enrich_assays and 'assay_confidence_score' in df.columns:
            logger.info(f"{uniprot_id}: Confidence scores - mean: {df['assay_confidence_score'].mean():.1f}, coverage: {df['assay_confidence_score'].notna().sum()}/{len(df)}")
        
        return df
    
    except Exception as e:
        logging.error(f"Error querying {uniprot_id}: {e}")
        return None


def process_system_group(args):
    """
    Process a group of systems with the same UniProt ID.
    Worker function for parallel processing.
    """
    uniprot_id, systems_df, output_dir, standard_types, max_activities, enrich_assays, enrich_targets, enrich_molecules, counter, total = args
    
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
        
        # Query ChEMBL with enrichment
        logger.info(f"[{counter.value}/{total}] Querying ChEMBL for {uniprot_id} ({len(systems_df)} systems)")
        
        activities_df = query_chembl_activities_enriched(
            uniprot_id,
            standard_types=standard_types,
            max_activities=max_activities,
            enrich_assays=enrich_assays,
            enrich_targets=enrich_targets,
            enrich_molecules=enrich_molecules
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
            'standard_types': activities_df['standard_type'].value_counts().to_dict() if 'standard_type' in activities_df.columns else {},
            'total_columns': len(activities_df.columns),
            'enrichment': {
                'assays': enrich_assays,
                'targets': enrich_targets,
                'molecules': enrich_molecules
            },
            'date_retrieved': pd.Timestamp.now().isoformat()
        }
        
        # Add enrichment-specific stats
        if enrich_assays and 'assay_confidence_score' in activities_df.columns:
            metadata['assay_confidence_stats'] = {
                'mean': float(activities_df['assay_confidence_score'].mean()),
                'median': float(activities_df['assay_confidence_score'].median()),
                'coverage': int(activities_df['assay_confidence_score'].notna().sum())
            }
        
        metadata_file = output_path / f"{uniprot_id}_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        counter.value += 1
        logger.info(f"[{counter.value}/{total}] ✓ {uniprot_id}: {len(activities_df)} activities, {len(systems_df)} systems, {len(activities_df.columns)} columns")
        
        return {
            'uniprot_id': uniprot_id,
            'success': True,
            'num_systems': len(systems_df),
            'num_activities': len(activities_df),
            'num_compounds': int(activities_df['molecule_chembl_id'].nunique()),
            'num_assays': int(activities_df['assay_chembl_id'].nunique()),
            'num_columns': len(activities_df.columns)
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
    """Merge ChEMBL activities with PLINDER system information."""
    logger = logging.getLogger(__name__)
    logger.info("Merging ChEMBL activities with PLINDER systems...")
    
    df = pd.read_parquet(annotation_file)
    df_exploded = df.explode('uniprot_ids')
    df_exploded = df_exploded[df_exploded['uniprot_ids'].notna()]
    
    all_data = []
    
    for uniprot_id in df_exploded['uniprot_ids'].unique():
        activity_file = Path(output_dir) / f"uniprot_{uniprot_id}" / f"{uniprot_id}_chembl_activities.parquet"
        
        if not activity_file.exists():
            continue
        
        activities = pd.read_parquet(activity_file)
        systems = df_exploded[df_exploded['uniprot_ids'] == uniprot_id]['system_id'].unique()
        
        activities['mapped_uniprot_id'] = uniprot_id
        activities['plinder_systems'] = [systems.tolist()] * len(activities)
        activities['num_plinder_systems'] = len(systems)
        
        all_data.append(activities)
    
    if not all_data:
        logger.error("No activities to merge!")
        return
    
    combined = pd.concat(all_data, ignore_index=True)
    
    output_path = Path(final_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_parquet(output_path, index=False)
    
    logger.info(f"✓ Merged data saved to {output_path}")
    logger.info(f"  Total activities: {len(combined)}")
    logger.info(f"  Unique compounds: {combined['molecule_chembl_id'].nunique()}")
    logger.info(f"  Unique UniProt IDs: {combined['mapped_uniprot_id'].nunique()}")
    logger.info(f"  Total columns: {len(combined.columns)}")


def main():
    parser = argparse.ArgumentParser(
        description='Enriched parallel ChEMBL data retrieval for PLINDER systems'
    )
    parser.add_argument('--input', default='data/annotation_table_with_uniprot.parquet')
    parser.add_argument('--output', default='data/chembl_parallel_enriched')
    parser.add_argument('--final-output', default='data/chembl_activities_enriched.parquet')
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--chunk-size', type=int, default=None)
    parser.add_argument('--max-activities', type=int, default=3000)
    parser.add_argument('--standard-types', nargs='+', default=['IC50', 'Ki', 'Kd', 'EC50'])
    parser.add_argument('--enrich-assays', action='store_true', default=True, help='Enrich with assay details (RECOMMENDED)')
    parser.add_argument('--no-enrich-assays', dest='enrich_assays', action='store_false')
    parser.add_argument('--enrich-targets', action='store_true', help='Enrich with target details')
    parser.add_argument('--enrich-molecules', action='store_true', help='Enrich with molecule properties')
    parser.add_argument('--skip-merge', action='store_true')
    
    args = parser.parse_args()
    
    logger = setup_logging(args.output)
    
    logger.info("="*80)
    logger.info("Enriched Parallel ChEMBL Data Retrieval for PLINDER Systems")
    logger.info("="*80)
    logger.info(f"Enrichment options:")
    logger.info(f"  Assays: {args.enrich_assays} (RECOMMENDED)")
    logger.info(f"  Targets: {args.enrich_targets}")
    logger.info(f"  Molecules: {args.enrich_molecules}")
    
    num_workers = args.workers if args.workers else max(1, cpu_count() - 1)
    logger.info(f"Using {num_workers} parallel workers")
    
    logger.info(f"Loading annotation table from {args.input}")
    df = pd.read_parquet(args.input)
    logger.info(f"Loaded {len(df)} systems")
    
    df_exploded = df.explode('uniprot_ids')
    df_exploded = df_exploded[df_exploded['uniprot_ids'].notna()]
    
    grouped = df_exploded.groupby('uniprot_ids')
    unique_uniprot_ids = list(grouped.groups.keys())
    logger.info(f"Found {len(unique_uniprot_ids)} unique UniProt IDs")
    
    work_items = []
    for uniprot_id, group_df in grouped:
        work_items.append((
            uniprot_id,
            group_df.reset_index(drop=True),
            args.output,
            args.standard_types,
            args.max_activities,
            args.enrich_assays,
            args.enrich_targets,
            args.enrich_molecules
        ))
    
    if args.chunk_size:
        chunks = [work_items[i:i+args.chunk_size] for i in range(0, len(work_items), args.chunk_size)]
        logger.info(f"Processing in {len(chunks)} chunks of {args.chunk_size}")
    else:
        chunks = [work_items]
    
    all_results = []
    
    for chunk_idx, chunk in enumerate(chunks, 1):
        if len(chunks) > 1:
            logger.info(f"\nProcessing chunk {chunk_idx}/{len(chunks)} ({len(chunk)} UniProt IDs)")
        
        manager = Manager()
        counter = manager.Value('i', 0)
        total = len(chunk)
        
        chunk_with_counter = [item + (counter, total) for item in chunk]
        
        with Pool(processes=num_workers) as pool:
            chunk_results = pool.map(process_system_group, chunk_with_counter)
        
        all_results.extend(chunk_results)
        
        if chunk_idx < len(chunks):
            time.sleep(5)
    
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
        total_columns = np.mean([r.get('num_columns', 0) for r in all_results if r['success'] and 'num_columns' in r])
        logger.info(f"\nTotal activities retrieved: {total_activities}")
        logger.info(f"Average columns per dataset: {total_columns:.0f}")
    
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
            'standard_types': args.standard_types,
            'enrich_assays': args.enrich_assays,
            'enrich_targets': args.enrich_targets,
            'enrich_molecules': args.enrich_molecules
        }
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"\n✓ Summary saved to {summary_file}")
    
    if not args.skip_merge:
        logger.info("\n" + "="*80)
        logger.info("MERGING RESULTS")
        logger.info("="*80)
        merge_activities_with_systems(args.input, args.output, args.final_output)
    
    logger.info("\n✓ All done!")


if __name__ == '__main__':
    main()
