#!/usr/bin/env python3
"""
Test script for enriched ChEMBL data retrieval using multiple clients.

This script demonstrates the correct approach:
1. Use activity client to get bioactivity data
2. Use assay client to get detailed assay information
3. Use target client for target details
4. Use molecule client for compound properties
5. Merge all information into a comprehensive dataset
"""

import pandas as pd
from chembl_webresource_client.new_client import new_client
from collections import defaultdict
import time

# Initialize clients
activity_client = new_client.activity
assay_client = new_client.assay
target_client = new_client.target
molecule_client = new_client.molecule

# Test with EGFR (P00533)
TEST_UNIPROT_ID = "P00533"

print("="*80)
print("Enriched ChEMBL Data Retrieval Test")
print("="*80)
print(f"\nTesting with UniProt ID: {TEST_UNIPROT_ID} (EGFR)")

# Step 1: Get activities
print("\nStep 1: Retrieving activities...")
activities = activity_client.filter(
    target_organism='Homo sapiens',
    standard_type__in=['IC50', 'Ki', 'Kd', 'EC50'],
    target_components__accession=TEST_UNIPROT_ID
)[:100]  # Limit to 100 for testing

print(f"Retrieved {len(activities)} activities")

# Convert to dataframe
activities_list = []
for act in activities:
    activities_list.append(dict(act))

activities_df = pd.DataFrame(activities_list)
print(f"Activities dataframe shape: {activities_df.shape}")
print(f"Columns from activity client: {len(activities_df.columns)}")

# Step 2: Get unique assay IDs and retrieve detailed assay info
print("\nStep 2: Retrieving detailed assay information...")
unique_assay_ids = activities_df['assay_chembl_id'].unique()
print(f"Found {len(unique_assay_ids)} unique assays")

assay_details_list = []
failed_assays = []

for i, assay_id in enumerate(unique_assay_ids):
    try:
        assay = assay_client.get(assay_id)
        assay_dict = {
            'assay_chembl_id': assay_id,
            # Assay identification
            'assay_description': assay.get('description'),
            'assay_type': assay.get('assay_type'),
            'assay_type_description': assay.get('assay_type_description'),
            
            # Assay conditions
            'assay_organism': assay.get('assay_organism'),
            'assay_strain': assay.get('assay_strain'),
            'assay_tissue': assay.get('assay_tissue'),
            'assay_cell_type': assay.get('assay_cell_type'),
            'assay_subcellular_fraction': assay.get('assay_subcellular_fraction'),
            'assay_tax_id': assay.get('assay_tax_id'),
            'assay_test_type': assay.get('assay_test_type'),
            
            # Quality metrics
            'confidence_score': assay.get('confidence_score'),
            'confidence_description': assay.get('confidence_description'),
            'relationship_type': assay.get('relationship_type'),
            'relationship_description': assay.get('relationship_description'),
            
            # BAO ontology
            'bao_format': assay.get('bao_format'),
            'bao_label': assay.get('bao_label'),
            
            # Links to other entities
            'cell_chembl_id': assay.get('cell_chembl_id'),
            'tissue_chembl_id': assay.get('tissue_chembl_id'),
            'target_chembl_id': assay.get('target_chembl_id'),
            'document_chembl_id': assay.get('document_chembl_id'),
            
            # Source information
            'src_id': assay.get('src_id'),
            'src_assay_id': assay.get('src_assay_id'),
            
            # Additional details
            'assay_category': assay.get('assay_category'),
            'assay_parameters': str(assay.get('assay_parameters')) if assay.get('assay_parameters') else None,
            'assay_classifications': str(assay.get('assay_classifications')) if assay.get('assay_classifications') else None,
            'variant_sequence': assay.get('variant_sequence')
        }
        assay_details_list.append(assay_dict)
        
        if (i + 1) % 10 == 0:
            print(f"  Processed {i + 1}/{len(unique_assay_ids)} assays...")
            time.sleep(0.5)  # Rate limiting
    
    except Exception as e:
        print(f"  Error retrieving {assay_id}: {e}")
        failed_assays.append(assay_id)

assay_details_df = pd.DataFrame(assay_details_list)
print(f"\nAssay details dataframe shape: {assay_details_df.shape}")
print(f"Successfully retrieved: {len(assay_details_list)}/{len(unique_assay_ids)} assays")
if failed_assays:
    print(f"Failed: {len(failed_assays)} assays")

# Step 3: Get unique target IDs and retrieve target info
print("\nStep 3: Retrieving target information...")
unique_target_ids = activities_df['target_chembl_id'].dropna().unique()
print(f"Found {len(unique_target_ids)} unique targets")

target_details_list = []
for target_id in unique_target_ids[:10]:  # Limit for testing
    try:
        target = target_client.get(target_id)
        target_dict = {
            'target_chembl_id': target_id,
            'target_pref_name': target.get('pref_name'),
            'target_organism': target.get('organism'),
            'target_type': target.get('target_type'),
            'target_tax_id': target.get('tax_id'),
            'species_group_flag': target.get('species_group_flag'),
            # Component information
            'target_components': str(target.get('target_components')) if target.get('target_components') else None
        }
        target_details_list.append(target_dict)
    except Exception as e:
        print(f"  Error retrieving target {target_id}: {e}")

target_details_df = pd.DataFrame(target_details_list)
print(f"Target details dataframe shape: {target_details_df.shape}")

# Step 4: Get unique molecule IDs and retrieve compound info
print("\nStep 4: Retrieving molecule information...")
unique_molecule_ids = activities_df['molecule_chembl_id'].dropna().unique()
print(f"Found {len(unique_molecule_ids)} unique molecules")

molecule_details_list = []
for mol_id in unique_molecule_ids[:50]:  # Limit for testing
    try:
        mol = molecule_client.get(mol_id)
        mol_dict = {
            'molecule_chembl_id': mol_id,
            'molecule_pref_name': mol.get('pref_name'),
            'molecule_max_phase': mol.get('max_phase'),
            'molecule_type': mol.get('molecule_type'),
            'first_approval': mol.get('first_approval'),
            'therapeutic_flag': mol.get('therapeutic_flag'),
            'natural_product': mol.get('natural_product'),
            # Molecular properties
            'molecule_properties': str(mol.get('molecule_properties')) if mol.get('molecule_properties') else None,
            'molecule_structures': str(mol.get('molecule_structures')) if mol.get('molecule_structures') else None
        }
        molecule_details_list.append(mol_dict)
    except Exception as e:
        print(f"  Error retrieving molecule {mol_id}: {e}")

molecule_details_df = pd.DataFrame(molecule_details_list)
print(f"Molecule details dataframe shape: {molecule_details_df.shape}")

# Step 5: Merge all data
print("\nStep 5: Merging all data...")

# Merge with assay details
enriched_df = activities_df.merge(
    assay_details_df,
    on='assay_chembl_id',
    how='left',
    suffixes=('', '_assay_enriched')
)
print(f"After assay merge: {enriched_df.shape}")

# Merge with target details
enriched_df = enriched_df.merge(
    target_details_df,
    on='target_chembl_id',
    how='left',
    suffixes=('', '_target_enriched')
)
print(f"After target merge: {enriched_df.shape}")

# Merge with molecule details
enriched_df = enriched_df.merge(
    molecule_details_df,
    on='molecule_chembl_id',
    how='left',
    suffixes=('', '_molecule_enriched')
)
print(f"After molecule merge: {enriched_df.shape}")

# Step 6: Save results
output_file = 'data/chembl_API_test/enriched_activities_test.parquet'
enriched_df.to_parquet(output_file, index=False)
print(f"\n✓ Enriched data saved to {output_file}")

# Show comparison
print("\n" + "="*80)
print("DATA ENRICHMENT COMPARISON")
print("="*80)

print("\n1. Assay Information Enrichment:")
original_assay_cols = [c for c in activities_df.columns if 'assay' in c.lower()]
enriched_assay_cols = [c for c in enriched_df.columns if 'assay' in c.lower()]
print(f"   Original assay-related columns: {len(original_assay_cols)}")
print(f"   Enriched assay-related columns: {len(enriched_assay_cols)}")
print(f"   New columns added: {len(enriched_assay_cols) - len(original_assay_cols)}")

print("\n2. Key Assay Fields Coverage:")
key_fields = ['assay_organism', 'assay_tissue', 'assay_cell_type', 
              'confidence_score', 'bao_format', 'bao_label']
for field in key_fields:
    if field in enriched_df.columns:
        non_null = enriched_df[field].notna().sum()
        print(f"   {field}: {non_null}/{len(enriched_df)} ({100*non_null/len(enriched_df):.1f}%)")

print("\n3. Total Columns:")
print(f"   Original (activity only): {len(activities_df.columns)}")
print(f"   Enriched (all clients): {len(enriched_df.columns)}")
print(f"   Additional columns: {len(enriched_df.columns) - len(activities_df.columns)}")

print("\n" + "="*80)
print("SUCCESS! Enriched dataset ready with complete assay/target/molecule information")
print("="*80)

# Show sample of enriched data
print("\nSample of enriched data (first 3 rows, key columns):")
sample_cols = ['molecule_chembl_id', 'standard_type', 'standard_value', 
               'assay_chembl_id', 'assay_type', 'assay_organism', 
               'assay_cell_type', 'confidence_score', 'bao_label',
               'target_chembl_id', 'target_type']
available_cols = [c for c in sample_cols if c in enriched_df.columns]
print(enriched_df[available_cols].head(3).to_string())
