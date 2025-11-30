#!/usr/bin/env python3
"""
Parallel PDB Chain to UniProt ID mapping for PLINDER systems.

This script efficiently maps PDB chains to UniProt IDs using multiprocessing.
Optimized for processing large numbers of PLINDER system IDs.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time
import json
import requests
from multiprocessing import Pool, Manager, cpu_count
import argparse
import logging
from datetime import datetime
import sys
import signal


def setup_logging(output_dir):
    """Setup logging to both file and console."""
    log_dir = Path(output_dir) if output_dir else Path('.')
    log_dir.mkdir(parents=True, exist_ok=True)
    
    log_file = log_dir / f"uniprot_mapping_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(processName)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    return logging.getLogger(__name__)


def parse_plinder_system_id(system_id):
    """
    Parse PLINDER system ID into components.
    Format: <PDB ID>__<assembly>__<receptor chain>__<ligand chain>
    
    Example: 7eek__1__1.A__1.I
    """
    parts = system_id.split('__')
    
    if len(parts) == 4:
        pdb_id, assembly, receptor_chain, ligand_chain = parts
        # Extract just the chain letter (after the dot)
        receptor_chain_only = receptor_chain.split('.')[-1] if '.' in receptor_chain else receptor_chain
        
        return {
            'pdb_id': pdb_id.lower(),
            'receptor_chain': receptor_chain_only,
            'assembly': assembly,
            'ligand_chain': ligand_chain
        }
    else:
        # Fallback for non-standard format
        return {
            'pdb_id': system_id.split('__')[0].lower() if '__' in system_id else system_id.lower(),
            'receptor_chain': None,
            'assembly': None,
            'ligand_chain': None
        }


def map_pdb_chain_to_uniprot_api(pdb_id, chain_id):
    """
    Map PDB ID + specific chain to UniProt accessions using RCSB GraphQL API.
    
    This uses the GraphQL API which provides chain-specific mappings.
    """
    logger = logging.getLogger(__name__)
    
    if not chain_id:
        logger.warning(f"No chain_id provided for {pdb_id}")
        return []
    
    # Use GraphQL API to get chain to UniProt mappings
    query = '''
    {
      entry(entry_id: "%s") {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
            auth_asym_ids
            uniprot_ids
          }
        }
      }
    }
    ''' % pdb_id.upper()
    
    try:
        url = 'https://data.rcsb.org/graphql'
        response = requests.post(url, json={'query': query}, timeout=10)
        
        if response.status_code != 200:
            logger.warning(f"API returned status {response.status_code} for {pdb_id}")
            return []
        
        data = response.json()
        
        # Check for GraphQL errors
        if 'errors' in data:
            logger.warning(f"GraphQL errors for {pdb_id}: {data['errors']}")
            return []
            
        if 'data' in data and data['data'] and 'entry' in data['data']:
            entry = data['data']['entry']
            if entry and 'polymer_entities' in entry:
                # Find the entity that contains our chain
                for entity in entry['polymer_entities']:
                    if 'rcsb_polymer_entity_container_identifiers' in entity:
                        identifiers = entity['rcsb_polymer_entity_container_identifiers']
                        auth_asym_ids = identifiers.get('auth_asym_ids', [])
                        uniprot_ids = identifiers.get('uniprot_ids', [])
                        
                        # Check if this entity contains our chain
                        if chain_id.upper() in [c.upper() for c in auth_asym_ids]:
                            if uniprot_ids:
                                return list(set(uniprot_ids))
                
                logger.debug(f"Chain {chain_id} not found in {pdb_id}")
            else:
                logger.warning(f"No polymer_entities for {pdb_id}")
        else:
            logger.warning(f"Unexpected API response structure for {pdb_id}")
            
    except requests.exceptions.Timeout:
        logger.error(f"Timeout for {pdb_id} chain {chain_id}")
    except requests.exceptions.RequestException as e:
        logger.error(f"Request error for {pdb_id} chain {chain_id}: {e}")
    except Exception as e:
        logger.error(f"Unexpected error for {pdb_id} chain {chain_id}: {e}")
    
    return []


def load_checkpoint(checkpoint_file):
    """Load existing mapping results from checkpoint file."""
    if Path(checkpoint_file).exists():
        try:
            with open(checkpoint_file, 'r') as f:
                return json.load(f)
        except Exception as e:
            logging.getLogger(__name__).warning(f"Failed to load checkpoint: {e}")
    return {}


def save_checkpoint(results, checkpoint_file):
    """Save mapping results to checkpoint file."""
    try:
        Path(checkpoint_file).parent.mkdir(parents=True, exist_ok=True)
        with open(checkpoint_file, 'w') as f:
            json.dump(results, f, indent=2)
    except Exception as e:
        logging.getLogger(__name__).error(f"Failed to save checkpoint: {e}")


def process_pdb_chain_batch(args):
    """
    Process a batch of PDB + Chain combinations.
    
    This is the worker function for parallel processing.
    
    Parameters:
    -----------
    args : tuple
        (batch_data, progress_dict, total_batches, batch_idx, checkpoint_file)
    
    Returns:
    --------
    dict : Mapping results for this batch
    """
    batch_data, progress_dict, total_batches, batch_idx, checkpoint_file = args
    logger = logging.getLogger(__name__)
    
    results = {}
    successful = 0
    failed = 0
    
    for pdb_id, chain_id in batch_data:
        pdb_chain_key = f"{pdb_id}_{chain_id}"
        
        try:
            uniprot_ids = map_pdb_chain_to_uniprot_api(pdb_id, chain_id)
            results[pdb_chain_key] = uniprot_ids
            
            if uniprot_ids:
                successful += 1
                logger.debug(f"✓ {pdb_id} chain {chain_id} → {uniprot_ids}")
            else:
                failed += 1
                logger.debug(f"✗ {pdb_id} chain {chain_id} → no mapping")
            
            # Rate limiting - be nice to RCSB API
            time.sleep(0.2)
            
        except Exception as e:
            logger.error(f"Error mapping {pdb_id} chain {chain_id}: {e}")
            results[pdb_chain_key] = []
            failed += 1
    
    # Update progress
    progress_dict[batch_idx] = len(batch_data)
    completed = sum(progress_dict.values())
    
    logger.info(f"Batch {batch_idx}/{total_batches} complete - Progress: {completed} mappings total ({successful} success, {failed} failed in this batch)")
    
    return results


def parallel_map_pdb_chains(pdb_chain_combos, num_workers=4, batch_size=10, checkpoint_file=None):
    """
    Map PDB + Chain combinations to UniProt IDs in parallel.
    
    Parameters:
    -----------
    pdb_chain_combos : list of tuples
        List of (pdb_id, chain_id) tuples
    num_workers : int
        Number of parallel workers
    batch_size : int
        Number of items per batch
    checkpoint_file : str
        Path to checkpoint file for resume capability
    
    Returns:
    --------
    dict : Mapping of pdb_chain_key -> list of UniProt IDs
    """
    logger = logging.getLogger(__name__)
    logger.info(f"Starting parallel mapping with {num_workers} workers")
    logger.info(f"Total PDB+Chain combinations: {len(pdb_chain_combos)}")
    
    # Load existing checkpoint if available
    all_results = {}
    if checkpoint_file:
        all_results = load_checkpoint(checkpoint_file)
        if all_results:
            logger.info(f"Loaded {len(all_results)} existing mappings from checkpoint")
            # Filter out already processed combinations
            existing_keys = set(all_results.keys())
            pdb_chain_combos = [
                combo for combo in pdb_chain_combos 
                if f"{combo[0]}_{combo[1]}" not in existing_keys
            ]
            logger.info(f"Remaining PDB+Chain combinations to process: {len(pdb_chain_combos)}")
            
            if not pdb_chain_combos:
                logger.info("All mappings already complete!")
                return all_results
    
    # Create batches
    batches = [pdb_chain_combos[i:i+batch_size] for i in range(0, len(pdb_chain_combos), batch_size)]
    logger.info(f"Created {len(batches)} batches of size ~{batch_size}")
    
    # Create shared progress tracker
    manager = Manager()
    progress_dict = manager.dict()
    
    # Prepare work items
    work_items = []
    for idx, batch in enumerate(batches):
        work_items.append((batch, progress_dict, len(batches), idx, checkpoint_file))
    
    # Setup signal handler for graceful shutdown (SLURM sends SIGTERM)
    # Use a container to allow signal handler to access all_results
    results_container = {'data': all_results}
    
    def signal_handler(signum, frame):
        logger.warning(f"Received signal {signum}! Saving checkpoint...")
        if checkpoint_file and results_container['data']:
            save_checkpoint(results_container['data'], checkpoint_file)
            logger.info(f"Checkpoint saved: {len(results_container['data'])} mappings")
        sys.exit(1)
    
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)  # Ctrl+C
    
    # Process in parallel with periodic checkpointing
    batch_count = 0
    try:
        with Pool(processes=num_workers) as pool:
            for i, batch_result in enumerate(pool.imap_unordered(process_pdb_chain_batch, work_items)):
                # Update results
                all_results.update(batch_result)
                results_container['data'] = all_results
                batch_count += 1
                
                logger.info(f"Collected batch {i+1}/{len(work_items)} - Total unique mappings so far: {len(all_results)}")
                
                # Save checkpoint every 10 batches or at the end
                if checkpoint_file and (batch_count % 10 == 0 or i == len(work_items) - 1):
                    save_checkpoint(all_results, checkpoint_file)
                    logger.info(f"✓ Checkpoint saved: {len(all_results)} mappings")
    
    except (KeyboardInterrupt, SystemExit) as e:
        logger.warning(f"Process interrupted: {e}! Saving checkpoint...")
        if checkpoint_file and all_results:
            save_checkpoint(all_results, checkpoint_file)
            logger.info(f"Checkpoint saved: {len(all_results)} mappings")
        raise
    
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if checkpoint_file and all_results:
            save_checkpoint(all_results, checkpoint_file)
            logger.info(f"Emergency checkpoint saved: {len(all_results)} mappings")
        raise
    
    logger.info(f"Mapping complete! Total mappings: {len(all_results)}")
    
    # Count successes
    success_count = sum(1 for uniprots in all_results.values() if uniprots)
    logger.info(f"Successful mappings: {success_count}/{len(all_results)}")
    
    return all_results


def main():
    parser = argparse.ArgumentParser(
        description='Parallel PDB Chain to UniProt ID mapping for PLINDER systems'
    )
    parser.add_argument(
        '--input',
        required=True,
        help='Input parquet file with system_id column'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output parquet file with UniProt mappings'
    )
    parser.add_argument(
        '--workers',
        type=int,
        default=None,
        help='Number of parallel workers (default: CPU count - 1)'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=10,
        help='Number of mappings per batch'
    )
    parser.add_argument(
        '--test',
        type=int,
        default=None,
        help='Test mode: only process first N systems'
    )
    parser.add_argument(
        '--log-dir',
        default='logs',
        help='Directory for log files'
    )
    parser.add_argument(
        '--checkpoint',
        default=None,
        help='Checkpoint file for resume capability (default: auto-generated in output directory)'
    )
    parser.add_argument(
        '--no-checkpoint',
        action='store_true',
        help='Disable checkpoint saving'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.log_dir)
    
    logger.info("="*80)
    logger.info("Parallel PDB Chain to UniProt Mapping")
    logger.info("="*80)
    
    # Determine number of workers
    num_workers = args.workers if args.workers else max(1, cpu_count() - 1)
    logger.info(f"Using {num_workers} parallel workers")
    
    # Test API connectivity
    logger.info("\nTesting RCSB API connectivity...")
    try:
        test_response = requests.get('https://data.rcsb.org/graphql', timeout=5)
        logger.info(f"✓ API accessible (status: {test_response.status_code})")
    except requests.exceptions.Timeout:
        logger.error("✗ API timeout - network issue or firewall blocking access")
    except requests.exceptions.ConnectionError as e:
        logger.error(f"✗ Cannot connect to API: {e}")
    except Exception as e:
        logger.error(f"✗ API test failed: {e}")
    
    # Setup checkpoint file
    checkpoint_file = None
    if not args.no_checkpoint:
        if args.checkpoint:
            checkpoint_file = args.checkpoint
        else:
            # Auto-generate checkpoint filename
            output_path = Path(args.output)
            checkpoint_file = output_path.parent / f".{output_path.stem}_checkpoint.json"
        logger.info(f"Checkpoint file: {checkpoint_file}")
    else:
        logger.info("Checkpoint disabled")
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    df = pd.read_parquet(args.input)
    logger.info(f"Loaded {len(df)} systems")
    
    # Test mode
    if args.test:
        df = df.head(args.test)
        logger.info(f"TEST MODE: Processing only {len(df)} systems")
    
    # Parse system IDs
    logger.info("Parsing PLINDER system IDs...")
    parsed_data = df['system_id'].apply(parse_plinder_system_id)
    
    df['pdb_id'] = parsed_data.apply(lambda x: x['pdb_id'])
    df['receptor_chain'] = parsed_data.apply(lambda x: x['receptor_chain'])
    df['ligand_chain'] = parsed_data.apply(lambda x: x['ligand_chain'])
    df['assembly_id'] = parsed_data.apply(lambda x: x['assembly'])
    
    # Get unique PDB + Chain combinations
    unique_combos = df[['pdb_id', 'receptor_chain']].drop_duplicates()
    unique_combos = unique_combos[unique_combos['receptor_chain'].notna()]
    
    pdb_chain_list = list(unique_combos.itertuples(index=False, name=None))
    
    logger.info(f"Unique PDB+Chain combinations: {len(pdb_chain_list)}")
    logger.info(f"Sample: {pdb_chain_list[:5]}")
    
    # Test single mapping before starting parallel processing
    if pdb_chain_list:
        logger.info("\nTesting single mapping...")
        test_pdb, test_chain = pdb_chain_list[0]
        test_result = map_pdb_chain_to_uniprot_api(test_pdb, test_chain)
        logger.info(f"Test mapping {test_pdb} chain {test_chain} → {test_result}")
        if not test_result:
            logger.warning("Test mapping returned no results - check API connectivity or data")
    
    # Parallel mapping
    logger.info("\n" + "="*80)
    logger.info("Starting parallel mapping...")
    logger.info("="*80 + "\n")
    
    start_time = time.time()
    
    mapping_results = parallel_map_pdb_chains(
        pdb_chain_list,
        num_workers=num_workers,
        batch_size=args.batch_size,
        checkpoint_file=checkpoint_file
    )
    
    elapsed_time = time.time() - start_time
    
    logger.info(f"\nMapping completed in {elapsed_time:.2f} seconds")
    if len(pdb_chain_list) > 0:
        logger.info(f"Average time per mapping: {elapsed_time/len(pdb_chain_list):.2f} seconds")
    
    # Debug: Check what we got
    logger.info(f"\nDEBUG: mapping_results type: {type(mapping_results)}")
    logger.info(f"DEBUG: mapping_results length: {len(mapping_results)}")
    if mapping_results:
        logger.info(f"DEBUG: First 3 mappings: {dict(list(mapping_results.items())[:3])}")
    else:
        logger.error("WARNING: mapping_results is empty!")
    
    # Add results to dataframe
    df['pdb_chain_key'] = df['pdb_id'] + '_' + df['receptor_chain']
    df['uniprot_ids'] = df['pdb_chain_key'].map(mapping_results)
    df['uniprot_ids_str'] = df['uniprot_ids'].apply(lambda x: ', '.join(x) if x else '')
    
    # Summary
    logger.info("\n" + "="*80)
    logger.info("SUMMARY")
    logger.info("="*80)
    logger.info(f"Total systems: {len(df)}")
    logger.info(f"Systems with UniProt IDs: {df['uniprot_ids_str'].astype(bool).sum()}")
    logger.info(f"Systems without UniProt IDs: {(~df['uniprot_ids_str'].astype(bool)).sum()}")
    logger.info(f"Unique UniProt IDs found: {len(set([uid for uids in df['uniprot_ids'].dropna() for uid in uids if uids]))}")
    
    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path)
    logger.info(f"\n✓ Results saved to: {output_path}")
    
    # Save mapping as JSON for reference
    mapping_file = output_path.parent / f"{output_path.stem}_mapping.json"
    with open(mapping_file, 'w') as f:
        json.dump(mapping_results, f, indent=2)
    logger.info(f"✓ Mapping dictionary saved to: {mapping_file}")
    
    # Save summary
    summary = {
        'total_systems': len(df),
        'unique_pdb_chain_combos': len(pdb_chain_list),
        'systems_with_uniprot': int(df['uniprot_ids_str'].astype(bool).sum()),
        'systems_without_uniprot': int((~df['uniprot_ids_str'].astype(bool)).sum()),
        'unique_uniprot_ids': len(set([uid for uids in df['uniprot_ids'].dropna() for uid in uids if uids])),
        'processing_time_seconds': elapsed_time,
        'workers_used': num_workers,
        'batch_size': args.batch_size,
        'date': datetime.now().isoformat()
    }
    
    summary_file = output_path.parent / f"{output_path.stem}_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
    logger.info(f"✓ Summary saved to: {summary_file}")
    
    # Clean up checkpoint file if successful
    if checkpoint_file and Path(checkpoint_file).exists():
        try:
            Path(checkpoint_file).unlink()
            logger.info(f"✓ Checkpoint file removed (job complete)")
        except Exception as e:
            logger.warning(f"Could not remove checkpoint file: {e}")
    
    logger.info("\n✓ All done!")


if __name__ == '__main__':
    main()
