
import pandas as pd
from rdkit import Chem
import sys
import os

# Add DeepCoy to path to get dataset_info
DEEPCOY_PATH = os.path.join(os.getcwd(), 'external', 'DeepCoy')
sys.path.insert(0, DEEPCOY_PATH)
from utils import dataset_info

def check_molecules(parquet_path):
    print(f"Checking {parquet_path}")
    df = pd.read_parquet(parquet_path)
    smiles_list = df['canonical_smiles'].unique()
    
    dataset = 'zinc'
    dataset_atoms = dataset_info(dataset)['atom_types']
    supported_atoms = set()
    for atom_type in dataset_atoms:
        symbol = ''.join([c for c in atom_type if c.isalpha()])
        supported_atoms.add(symbol)
    
    print(f"Supported atoms in {dataset}: {sorted(list(supported_atoms))}")
    
    for i, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Invalid SMILES: {smiles}")
            continue
            
        num_heavy = mol.GetNumHeavyAtoms()
        atoms = set([a.GetSymbol() for a in mol.GetAtoms()])
        unsupported = atoms - supported_atoms
        
        status = "PASS"
        reasons = []
        if num_heavy < 10:
            status = "FAIL"
            reasons.append(f"Too small ({num_heavy} < 10)")
        if num_heavy > 80:
            status = "FAIL"
            reasons.append(f"Too large ({num_heavy} > 80)")
        if unsupported:
            status = "FAIL"
            reasons.append(f"Unsupported atoms: {unsupported}")
            
        print(f"Mol {i+1}: {smiles[:30]}... | Heavy: {num_heavy} | Atoms: {atoms} | {status} {reasons}")

if __name__ == "__main__":
    check_molecules('/project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity/uniprot_O43598/O43598_chembl_activities_filtered.parquet')
