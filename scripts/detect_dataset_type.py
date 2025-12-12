
import pandas as pd
from rdkit import Chem
import sys
import argparse

def detect_dataset(parquet_path):
    try:
        df = pd.read_parquet(parquet_path)
        # Check a sample of molecules for Phosphorus
        # If any molecule has P, we should use the phosphorus model
        # We check all unique SMILES to be safe, but could limit if too slow
        unique_smiles = df['canonical_smiles'].dropna().unique()
        
        has_phosphorus = False
        for smiles in unique_smiles:
            if 'P' in smiles or 'p' in smiles: # Quick string check first
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
                    if 'P' in atoms:
                        has_phosphorus = True
                        break
        
        if has_phosphorus:
            print("zinc_phosphorus")
        else:
            print("zinc")
            
    except Exception as e:
        # Fallback to default if error
        sys.stderr.write(f"Error detecting dataset: {e}\n")
        print("zinc")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input_parquet', help='Path to input parquet file')
    args = parser.parse_args()
    detect_dataset(args.input_parquet)
