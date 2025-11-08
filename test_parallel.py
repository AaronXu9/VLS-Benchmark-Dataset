#!/usr/bin/env python3
"""
Test script for parallel ChEMBL retrieval - local testing before CARC deployment.

This tests the parallel retrieval with a small subset to verify everything works.
"""

import pandas as pd
from pathlib import Path
import subprocess
import sys

def main():
    print("="*70)
    print("Testing Parallel ChEMBL Retrieval (Local)")
    print("="*70)
    
    # Check input file exists
    input_file = Path('data/annotation_table_with_uniprot.parquet')
    if not input_file.exists():
        print(f"❌ Input file not found: {input_file}")
        print("   Run extract_uniprot_ids.ipynb first!")
        sys.exit(1)
    
    print(f"✓ Input file found: {input_file}")
    
    # Load and check data
    df = pd.read_parquet(input_file)
    print(f"✓ Loaded {len(df)} systems")
    
    # Check for UniProt IDs
    has_uniprot = df['uniprot_ids'].notna().sum()
    print(f"✓ Systems with UniProt IDs: {has_uniprot}/{len(df)}")
    
    if has_uniprot == 0:
        print("❌ No UniProt IDs found!")
        sys.exit(1)
    
    # Create test subset
    test_df = df[df['uniprot_ids'].notna()].head(20)
    test_file = Path('data/annotation_table_test.parquet')
    test_df.to_parquet(test_file)
    print(f"✓ Created test subset: {len(test_df)} systems")
    
    # Get unique UniProt IDs in test
    df_exploded = test_df.explode('uniprot_ids')
    unique_uniprots = df_exploded['uniprot_ids'].nunique()
    print(f"✓ Test contains {unique_uniprots} unique UniProt IDs")
    
    # Run test
    print("\n" + "="*70)
    print("Running parallel retrieval test (4 workers)...")
    print("="*70 + "\n")
    
    cmd = [
        sys.executable,
        'retrieve_chembl_parallel.py',
        '--input', str(test_file),
        '--output', 'data/chembl_test',
        '--final-output', 'data/chembl_test_merged.parquet',
        '--workers', '4',
        '--max-activities', '100'  # Limit for testing
    ]
    
    result = subprocess.run(cmd)
    
    if result.returncode == 0:
        print("\n" + "="*70)
        print("✓ Test completed successfully!")
        print("="*70)
        
        # Check outputs
        summary_file = Path('data/chembl_test/retrieval_summary.json')
        merged_file = Path('data/chembl_test_merged.parquet')
        
        if summary_file.exists():
            print(f"✓ Summary file created: {summary_file}")
            
            import json
            with open(summary_file) as f:
                summary = json.load(f)
            
            print(f"\n  Results:")
            print(f"    - Total UniProt IDs: {summary['total_uniprot_ids']}")
            print(f"    - Successful: {summary['successful']}")
            print(f"    - Failed: {summary['failed']}")
        
        if merged_file.exists():
            print(f"✓ Merged file created: {merged_file}")
            
            merged_df = pd.read_parquet(merged_file)
            print(f"\n  Merged data:")
            print(f"    - Total activities: {len(merged_df)}")
            print(f"    - Unique compounds: {merged_df['molecule_chembl_id'].nunique()}")
            print(f"    - Columns: {len(merged_df.columns)}")
        
        print("\n" + "="*70)
        print("Ready for CARC deployment!")
        print("="*70)
        print("\nNext steps:")
        print("  1. Review test outputs in data/chembl_test/")
        print("  2. If satisfied, transfer files to CARC")
        print("  3. Submit job with: sbatch run_chembl_parallel.slurm")
        
    else:
        print("\n❌ Test failed! Check errors above.")
        sys.exit(1)


if __name__ == '__main__':
    main()
