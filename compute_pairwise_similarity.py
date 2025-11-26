#!/usr/bin/env python3
"""
Compute pairwise PLI similarity scores for systems in matched_annotation_table.parquet.

This script computes similarity scores between all pairs of protein-ligand systems,
including:
1. Tanimoto (2D fingerprint similarity)
2. SuCOS Protein-Aligned (3D overlap after Foldseek protein superposition)
3. SuCOS Ligand-Aligned (3D overlap after optimal ligand alignment)

Usage:
    # Compute similarities for all pairs (very slow for large datasets!)
    python compute_pairwise_similarity.py --all-pairs
    
    # Compute similarities within the same PDB
    python compute_pairwise_similarity.py --same-pdb
    
    # Compute similarities for a specific query system against all others
    python compute_pairwise_similarity.py --query-system 4grb__1__1.A__1.C
    
    # Compute similarities between two specific systems
    python compute_pairwise_similarity.py --query-system 4grb__1__1.A__1.C --target-system 6gra__1__1.A__1.B
    
    # Use pre-computed Foldseek alignments (skip protein-aligned scoring if not available)
    python compute_pairwise_similarity.py --same-pdb --skip-protein-aligned
    
    # Parallel processing with multiple workers
    python compute_pairwise_similarity.py --same-pdb --n-workers 4
"""

import argparse
import logging
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
from typing import List, Optional, Tuple

import pandas as pd
from tqdm import tqdm 

from pli_similarity_scorer import (
    PLISimilarityScorer,
    load_annotation_table,
    SimilarityScore,
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
LOG = logging.getLogger(__name__)


def get_unique_system_ids(annotation_table: pd.DataFrame) -> List[str]:
    """Get unique system IDs from annotation table."""
    return annotation_table["system_id"].unique().tolist()


def get_system_pairs_same_pdb(annotation_table: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Generate pairs of systems from the same PDB entry.
    This is useful for finding similar binding modes within a structure.
    """
    pairs = []
    
    # Group by PDB ID (first 4 characters of system_id)
    annotation_table = annotation_table.copy()
    annotation_table["pdb_id"] = annotation_table["system_id"].str[:4]
    
    for pdb_id, group in annotation_table.groupby("pdb_id"):
        system_ids = group["system_id"].unique().tolist()
        if len(system_ids) >= 2:
            pairs.extend(combinations(system_ids, 2))
    
    return pairs


def get_system_pairs_all(annotation_table: pd.DataFrame) -> List[Tuple[str, str]]:
    """
    Generate all possible pairs of systems.
    WARNING: This can be extremely large! N*(N-1)/2 pairs.
    """
    system_ids = get_unique_system_ids(annotation_table)
    return list(combinations(system_ids, 2))


def get_system_pairs_query_vs_all(
    annotation_table: pd.DataFrame,
    query_system: str,
) -> List[Tuple[str, str]]:
    """Generate pairs of query system vs all other systems."""
    system_ids = get_unique_system_ids(annotation_table)
    return [(query_system, s) for s in system_ids if s != query_system]


def score_pair_worker(args: Tuple) -> dict:
    """Worker function for parallel scoring."""
    (
        system_1,
        system_2,
        systems_dir,
        foldseek_cache_dir,
        compute_protein_aligned,
    ) = args
    
    # Create a new scorer instance for this worker
    scorer = PLISimilarityScorer(
        systems_dir=systems_dir,
        foldseek_cache_dir=foldseek_cache_dir,
    )
    
    try:
        score = scorer.score_pair(
            system_1,
            system_2,
            compute_protein_aligned=compute_protein_aligned,
        )
        return {
            "query_system": score.query_system,
            "target_system": score.target_system,
            "tanimoto": score.tanimoto,
            "sucos_protein_aligned": score.sucos_protein_aligned,
            "sucos_ligand_aligned": score.sucos_ligand_aligned,
            "shape_similarity": score.shape_similarity,
            "alignment_rmsd": score.alignment_rmsd,
            "foldseek_lddt": score.foldseek_lddt,
            "error": None,
        }
    except Exception as e:
        return {
            "query_system": system_1,
            "target_system": system_2,
            "error": str(e),
        }


def main():
    parser = argparse.ArgumentParser(
        description="Compute pairwise PLI similarity scores",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    
    # Input/output paths
    parser.add_argument(
        "--annotation-table",
        type=Path,
        default=Path("data/matched_annotation_table.parquet"),
        help="Path to matched_annotation_table.parquet",
    )
    parser.add_argument(
        "--systems-dir",
        type=Path,
        default=Path("/mnt/katritch_lab2/aoxu/2024-06/v2/systems"),
        help="Directory containing system folders",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("similarity_scores/pairwise_scores.parquet"),
        help="Output file for similarity scores",
    )
    parser.add_argument(
        "--foldseek-cache-dir",
        type=Path,
        default=Path("similarity_scores/foldseek_cache"),
        help="Directory for cached Foldseek alignments",
    )
    
    # Pair selection modes
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--all-pairs",
        action="store_true",
        help="Compute all pairwise similarities (WARNING: very slow!)",
    )
    mode_group.add_argument(
        "--same-pdb",
        action="store_true",
        help="Only compute similarities between systems from the same PDB",
    )
    mode_group.add_argument(
        "--query-system",
        type=str,
        help="Query system ID to compare against all others (or specific target)",
    )
    
    parser.add_argument(
        "--target-system",
        type=str,
        help="Specific target system ID (requires --query-system)",
    )
    
    # Scoring options
    parser.add_argument(
        "--skip-protein-aligned",
        action="store_true",
        help="Skip protein-aligned scoring (much faster, no Foldseek needed)",
    )
    
    # Performance options
    parser.add_argument(
        "--n-workers",
        type=int,
        default=1,
        help="Number of parallel workers (default: 1, sequential)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=1000,
        help="Save intermediate results every N pairs",
    )
    parser.add_argument(
        "--max-pairs",
        type=int,
        default=None,
        help="Maximum number of pairs to process (for testing)",
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.target_system and not args.query_system:
        parser.error("--target-system requires --query-system")
    
    # Load annotation table
    LOG.info(f"Loading annotation table from {args.annotation_table}")
    annotation_table = load_annotation_table(args.annotation_table)
    LOG.info(f"Loaded {len(annotation_table)} systems")
    
    # Determine pairs to score
    if args.all_pairs:
        pairs = get_system_pairs_all(annotation_table)
        LOG.info(f"Computing ALL pairwise similarities: {len(pairs)} pairs")
    elif args.same_pdb:
        pairs = get_system_pairs_same_pdb(annotation_table)
        LOG.info(f"Computing same-PDB similarities: {len(pairs)} pairs")
    elif args.query_system:
        if args.target_system:
            pairs = [(args.query_system, args.target_system)]
            LOG.info(f"Computing similarity between {args.query_system} and {args.target_system}")
        else:
            pairs = get_system_pairs_query_vs_all(annotation_table, args.query_system)
            LOG.info(f"Computing {args.query_system} vs all: {len(pairs)} pairs")
    
    # Limit pairs if requested
    if args.max_pairs and len(pairs) > args.max_pairs:
        pairs = pairs[:args.max_pairs]
        LOG.info(f"Limited to {args.max_pairs} pairs")
    
    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.foldseek_cache_dir.mkdir(parents=True, exist_ok=True)
    
    compute_protein_aligned = not args.skip_protein_aligned
    
    # Score pairs
    results = []
    
    if args.n_workers <= 1:
        # Sequential processing
        scorer = PLISimilarityScorer(
            systems_dir=args.systems_dir,
            annotation_table=annotation_table,
            foldseek_cache_dir=args.foldseek_cache_dir,
        )
        
        for i, (system_1, system_2) in enumerate(tqdm(pairs, desc="Scoring pairs")):
            try:
                score = scorer.score_pair(
                    system_1,
                    system_2,
                    compute_protein_aligned=compute_protein_aligned,
                )
                results.append({
                    "query_system": score.query_system,
                    "target_system": score.target_system,
                    "tanimoto": score.tanimoto,
                    "sucos_protein_aligned": score.sucos_protein_aligned,
                    "sucos_ligand_aligned": score.sucos_ligand_aligned,
                    "shape_similarity": score.shape_similarity,
                    "alignment_rmsd": score.alignment_rmsd,
                    "foldseek_lddt": score.foldseek_lddt,
                })
            except Exception as e:
                LOG.warning(f"Error scoring {system_1} vs {system_2}: {e}")
                results.append({
                    "query_system": system_1,
                    "target_system": system_2,
                    "error": str(e),
                })
            
            # Save intermediate results
            if (i + 1) % args.batch_size == 0:
                df = pd.DataFrame(results)
                df.to_parquet(args.output, index=False)
                LOG.info(f"Saved {len(results)} results to {args.output}")
    
    else:
        # Parallel processing
        worker_args = [
            (s1, s2, args.systems_dir, args.foldseek_cache_dir, compute_protein_aligned)
            for s1, s2 in pairs
        ]
        
        with ProcessPoolExecutor(max_workers=args.n_workers) as executor:
            futures = {executor.submit(score_pair_worker, arg): arg for arg in worker_args}
            
            for future in tqdm(as_completed(futures), total=len(futures), desc="Scoring pairs"):
                result = future.result()
                results.append(result)
                
                # Save intermediate results
                if len(results) % args.batch_size == 0:
                    df = pd.DataFrame(results)
                    df.to_parquet(args.output, index=False)
                    LOG.info(f"Saved {len(results)} results to {args.output}")
    
    # Save final results
    df = pd.DataFrame(results)
    df.to_parquet(args.output, index=False)
    LOG.info(f"Saved {len(results)} final results to {args.output}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total pairs scored: {len(df)}")
    
    if "error" in df.columns:
        errors = df["error"].notna().sum()
        print(f"Errors: {errors}")
    
    if "tanimoto" in df.columns:
        valid = df["tanimoto"].notna()
        print(f"\nTanimoto (2D):")
        print(f"  Valid scores: {valid.sum()}")
        print(f"  Mean: {df.loc[valid, 'tanimoto'].mean():.4f}")
        print(f"  Median: {df.loc[valid, 'tanimoto'].median():.4f}")
    
    if "sucos_ligand_aligned" in df.columns:
        valid = df["sucos_ligand_aligned"].notna()
        print(f"\nSuCOS Ligand-Aligned (3D):")
        print(f"  Valid scores: {valid.sum()}")
        print(f"  Mean: {df.loc[valid, 'sucos_ligand_aligned'].mean():.4f}")
        print(f"  Median: {df.loc[valid, 'sucos_ligand_aligned'].median():.4f}")
    
    if "sucos_protein_aligned" in df.columns and not args.skip_protein_aligned:
        valid = df["sucos_protein_aligned"].notna()
        print(f"\nSuCOS Protein-Aligned (3D):")
        print(f"  Valid scores: {valid.sum()}")
        if valid.sum() > 0:
            print(f"  Mean: {df.loc[valid, 'sucos_protein_aligned'].mean():.4f}")
            print(f"  Median: {df.loc[valid, 'sucos_protein_aligned'].median():.4f}")
    
    print(f"\nResults saved to: {args.output}")
    

if __name__ == "__main__":
    main()
