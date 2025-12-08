#!/usr/bin/env python
"""
Convert ChEMBL affinity data (parquet) to DeepCoy input format (JSON).

Usage:
    python prepare_deepcoy_input.py --input_parquet FILE --output_json FILE [options]
"""

import sys
import os
import argparse

# Add DeepCoy to path
DEEPCOY_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'external', 'DeepCoy')
sys.path.insert(0, DEEPCOY_PATH)

import json
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops

# Import DeepCoy utilities
from utils import bond_dict, dataset_info, to_graph_mol


def smiles_to_deepcoy_format(smiles, dataset='zinc'):
    """
    Convert a SMILES string to DeepCoy graph format.
    
    For generation (not training), we set smiles_in = smiles_out and 
    graph_in = graph_out since we want to generate decoys from the same molecule.
    
    Args:
        smiles: SMILES string
        dataset: 'zinc' or 'zinc_phosphorus'
    
    Returns:
        dict with DeepCoy format or None if conversion fails
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        
        # Convert to graph representation
        nodes, edges = to_graph_mol(mol, dataset)
        
        if len(edges) <= 0:
            return None
            
        # For generation, input and output are the same
        return {
            'graph_in': edges,
            'graph_out': edges,
            'node_features_in': nodes,
            'node_features_out': nodes,
            'smiles_out': canonical_smiles,
            'smiles_in': canonical_smiles,
            'v_to_keep': []  # Empty for generation
        }
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        return None


def filter_molecule(smiles, min_heavy_atoms=10, max_heavy_atoms=80, dataset='zinc'):
    """
    Check if molecule passes filters for DeepCoy.
    
    Args:
        smiles: SMILES string
        min_heavy_atoms: Minimum number of heavy atoms
        max_heavy_atoms: Maximum number of heavy atoms  
        dataset: Dataset type for atom type checking
    
    Returns:
        bool: True if molecule passes filters
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        # Check heavy atom count
        num_heavy = mol.GetNumHeavyAtoms()
        if num_heavy < min_heavy_atoms or num_heavy > max_heavy_atoms:
            return False
        
        # Check if molecule contains only supported atom types
        dataset_atoms = dataset_info(dataset)['atom_types']
        # Extract base atom symbols from dataset (e.g., 'C4(0)' -> 'C')
        supported_atoms = set()
        for atom_type in dataset_atoms:
            symbol = ''.join([c for c in atom_type if c.isalpha()])
            supported_atoms.add(symbol)
        
        for atom in mol.GetAtoms():
            if atom.GetSymbol() not in supported_atoms:
                return False
        
        return True
    except:
        return False


def convert_parquet_to_deepcoy(input_parquet, output_json, smiles_column='canonical_smiles',
                                min_heavy_atoms=10, max_heavy_atoms=80, dataset='zinc',
                                max_molecules=0):
    """
    Convert a parquet file with SMILES to DeepCoy JSON format.
    
    Args:
        input_parquet: Path to input parquet file
        output_json: Path to output JSON file
        smiles_column: Column name containing SMILES
        min_heavy_atoms: Minimum heavy atoms filter
        max_heavy_atoms: Maximum heavy atoms filter
        dataset: 'zinc' or 'zinc_phosphorus'
        max_molecules: Maximum number of molecules to process (0 = no limit)
    
    Returns:
        tuple: (num_processed, num_failed, num_filtered)
    """
    # Load parquet
    print(f"Loading {input_parquet}...")
    df = pd.read_parquet(input_parquet)
    
    if smiles_column not in df.columns:
        raise ValueError(f"Column '{smiles_column}' not found in parquet. Available columns: {df.columns.tolist()}")
    
    # Get unique SMILES
    unique_smiles = df[smiles_column].dropna().unique()
    print(f"Found {len(unique_smiles)} unique SMILES")
    
    # Limit molecules if max_molecules is set
    if max_molecules > 0 and len(unique_smiles) > max_molecules:
        print(f"Limiting to {max_molecules} molecules (test mode)")
        unique_smiles = unique_smiles[:max_molecules]
    
    processed_data = []
    num_filtered = 0
    num_failed = 0
    
    for i, smiles in enumerate(unique_smiles):
        if i % 100 == 0:
            print(f"Processing {i}/{len(unique_smiles)}...", end='\r')
        
        # Apply filters
        if not filter_molecule(smiles, min_heavy_atoms, max_heavy_atoms, dataset):
            num_filtered += 1
            continue
        
        # Convert to DeepCoy format
        mol_data = smiles_to_deepcoy_format(smiles, dataset)
        if mol_data is None:
            num_failed += 1
            continue
        
        processed_data.append(mol_data)
    
    print(f"\nProcessing complete!")
    print(f"  Total unique SMILES: {len(unique_smiles)}")
    print(f"  Successfully converted: {len(processed_data)}")
    print(f"  Filtered out: {num_filtered}")
    print(f"  Failed to convert: {num_failed}")
    
    # Save to JSON
    print(f"Saving to {output_json}...")
    os.makedirs(os.path.dirname(output_json) if os.path.dirname(output_json) else '.', exist_ok=True)
    with open(output_json, 'w') as f:
        json.dump(processed_data, f)
    
    print(f"Done! Saved {len(processed_data)} molecules.")
    return len(processed_data), num_failed, num_filtered


if __name__ == "__main__":
    import re
    
    parser = argparse.ArgumentParser(description='Convert ChEMBL parquet to DeepCoy JSON format')
    parser.add_argument('--input_parquet', required=True, help='Path to input parquet file')
    parser.add_argument('--output_json', help='Path to output JSON file (auto-generated if not provided)')
    parser.add_argument('--smiles_column', default='canonical_smiles', help='Column name containing SMILES')
    parser.add_argument('--uniprot_id', help='UniProt ID for naming (extracted from path if not provided)')
    parser.add_argument('--min_heavy_atoms', type=int, default=10, help='Minimum heavy atoms')
    parser.add_argument('--max_heavy_atoms', type=int, default=80, help='Maximum heavy atoms')
    parser.add_argument('--dataset', default='zinc', choices=['zinc', 'zinc_phosphorus'], help='Dataset type')
    parser.add_argument('--max_molecules', type=int, default=0, help='Maximum molecules to process (0 = no limit, for testing)')
    
    args = parser.parse_args()
    
    output_json = args.output_json
    if not output_json:
        # Auto-generate output name
        uniprot_id = args.uniprot_id
        if not uniprot_id:
            match = re.search(r'uniprot_([A-Z0-9]+)', args.input_parquet)
            if match:
                uniprot_id = match.group(1)
            else:
                uniprot_id = 'unknown'
        output_json = f"molecules_{uniprot_id}_actives.json"
    
    convert_parquet_to_deepcoy(
        input_parquet=args.input_parquet,
        output_json=output_json,
        smiles_column=args.smiles_column,
        min_heavy_atoms=args.min_heavy_atoms,
        max_heavy_atoms=args.max_heavy_atoms,
        dataset=args.dataset,
        max_molecules=args.max_molecules
    )
