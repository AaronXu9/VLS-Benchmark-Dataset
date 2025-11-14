#!/usr/bin/env python3
"""
Test script to explore ChEMBL API clients and their capabilities.

This script tests different ChEMBL clients to understand:
1. What information each client provides
2. How to properly retrieve assay-level information
3. The relationship between different data types (activities, assays, targets, molecules)
"""

import pandas as pd
from chembl_webresource_client.new_client import new_client
import json
from pprint import pprint

# Initialize different ChEMBL clients
# Test which clients are available
available_clients = {}
client_names = ['activity', 'assay', 'target', 'molecule', 'mechanism', 'drug', 
                'document', 'cell_line', 'tissue', 'protein_class', 'source', 
                'binding_site', 'target_component', 'chembl_id_lookup', 'image',
                'atc_class', 'drug_indication', 'go_slim']

print("Testing available clients...")
for client_name in client_names:
    try:
        client = getattr(new_client, client_name)
        available_clients[client_name] = client
        print(f"  ✓ {client_name}")
    except AttributeError:
        print(f"  ✗ {client_name} (not available)")

print(f"\nFound {len(available_clients)} available clients\n")

# Assign commonly used clients
activity_client = available_clients.get('activity')
assay_client = available_clients.get('assay')
target_client = available_clients.get('target')
molecule_client = available_clients.get('molecule')
mechanism_client = available_clients.get('mechanism')
document_client = available_clients.get('document')
target_component_client = available_clients.get('target_component')

print("="*80)
print("ChEMBL API Clients Exploration")
print("="*80)

# Test with a known UniProt ID (e.g., EGFR - P00533)
TEST_UNIPROT_ID = "P00533"  # EGFR
TEST_MOLECULE_CHEMBL_ID = "CHEMBL25"  # Erlotinib (EGFR inhibitor)

print(f"\nTest UniProt ID: {TEST_UNIPROT_ID} (EGFR)")
print(f"Test Molecule: {TEST_MOLECULE_CHEMBL_ID} (Erlotinib)")

# =============================================================================
# 1. ACTIVITY CLIENT - What we're currently using
# =============================================================================
print("\n" + "="*80)
print("1. ACTIVITY CLIENT")
print("="*80)

activities = activity_client.filter(
    target_organism='Homo sapiens',
    standard_type__in=['IC50', 'Ki', 'Kd', 'EC50'],
    target_components__accession=TEST_UNIPROT_ID
).only(['molecule_chembl_id', 'standard_type', 'standard_value', 
        'assay_chembl_id', 'assay_description', 'assay_type'])[:5]

print(f"\nFound {len(activities)} activities (limited to 5)")
if activities:
    print("\nFirst activity fields:")
    first_activity = activities[0]
    pprint(dict(first_activity))
    
    print("\nAll available fields in activity object:")
    pprint(sorted(first_activity.keys()))

# =============================================================================
# 2. ASSAY CLIENT - For detailed assay information
# =============================================================================
print("\n" + "="*80)
print("2. ASSAY CLIENT")
print("="*80)

# Get assay IDs from activities
if activities:
    test_assay_id = activities[0]['assay_chembl_id']
    print(f"\nTesting with assay: {test_assay_id}")
    
    # Get detailed assay information
    assay_details = assay_client.get(test_assay_id)
    
    print("\nDetailed assay information:")
    pprint(dict(assay_details))
    
    print("\nAll available fields in assay object:")
    pprint(sorted(assay_details.keys()))
    
    # Try filtering assays by target
    print(f"\nSearching assays for UniProt ID: {TEST_UNIPROT_ID}")
    assays = assay_client.filter(
        target_organism='Homo sapiens',
        target_components__accession=TEST_UNIPROT_ID
    )[:10]
    
    print(f"Found {len(assays)} assays (limited to 10)")
    
    # Compare assay info from activity vs assay client
    print("\n" + "-"*80)
    print("COMPARISON: Assay info from Activity vs Assay Client")
    print("-"*80)
    
    print("\nFrom Activity client:")
    print(f"  assay_description: {activities[0].get('assay_description')}")
    print(f"  assay_type: {activities[0].get('assay_type')}")
    print(f"  assay_organism: {activities[0].get('assay_organism')}")
    
    print("\nFrom Assay client (has more details):")
    print(f"  description: {assay_details.get('description')}")
    print(f"  assay_type: {assay_details.get('assay_type')}")
    print(f"  assay_organism: {assay_details.get('assay_organism')}")
    print(f"  assay_strain: {assay_details.get('assay_strain')}")
    print(f"  assay_tissue: {assay_details.get('assay_tissue')}")
    print(f"  assay_cell_type: {assay_details.get('assay_cell_type')}")
    print(f"  assay_subcellular_fraction: {assay_details.get('assay_subcellular_fraction')}")
    print(f"  confidence_score: {assay_details.get('confidence_score')}")
    print(f"  src_id: {assay_details.get('src_id')}")
    print(f"  src_assay_id: {assay_details.get('src_assay_id')}")
    print(f"  bao_format: {assay_details.get('bao_format')}")
    print(f"  bao_label: {assay_details.get('bao_label')}")
    print(f"  cell_chembl_id: {assay_details.get('cell_chembl_id')}")
    print(f"  tissue_chembl_id: {assay_details.get('tissue_chembl_id')}")
    print(f"  variant_sequence: {assay_details.get('variant_sequence')}")

# =============================================================================
# 3. TARGET CLIENT - For detailed target information
# =============================================================================
print("\n" + "="*80)
print("3. TARGET CLIENT")
print("="*80)

# Search targets by UniProt ID
targets = target_client.filter(
    target_components__accession=TEST_UNIPROT_ID
)

print(f"\nFound {len(targets)} targets")
if targets:
    target = targets[0]
    print("\nFirst target details:")
    pprint(dict(target))
    
    print("\nAll available fields in target object:")
    pprint(sorted(target.keys()))

# =============================================================================
# 4. MOLECULE CLIENT - For detailed compound information
# =============================================================================
print("\n" + "="*80)
print("4. MOLECULE CLIENT")
print("="*80)

molecule = molecule_client.get(TEST_MOLECULE_CHEMBL_ID)

print(f"\nMolecule details for {TEST_MOLECULE_CHEMBL_ID}:")
pprint(dict(molecule))

print("\nAll available fields in molecule object:")
pprint(sorted(molecule.keys()))

# =============================================================================
# 5. DOCUMENT CLIENT - For publication information
# =============================================================================
print("\n" + "="*80)
print("5. DOCUMENT CLIENT")
print("="*80)

if activities:
    doc_id = activities[0].get('document_chembl_id')
    if doc_id:
        print(f"\nTesting with document: {doc_id}")
        doc = document_client.get(doc_id)
        
        print("\nDocument details:")
        pprint(dict(doc))
        
        print("\nAll available fields in document object:")
        pprint(sorted(doc.keys()))

# =============================================================================
# 6. MECHANISM CLIENT - For drug mechanism information
# =============================================================================
print("\n" + "="*80)
print("6. MECHANISM CLIENT")
print("="*80)

mechanisms = mechanism_client.filter(molecule_chembl_id=TEST_MOLECULE_CHEMBL_ID)

print(f"\nFound {len(mechanisms)} mechanisms for {TEST_MOLECULE_CHEMBL_ID}")
if mechanisms:
    print("\nFirst mechanism:")
    pprint(dict(mechanisms[0]))
    
    print("\nAll available fields in mechanism object:")
    pprint(sorted(mechanisms[0].keys()))

# =============================================================================
# 7. TARGET COMPONENT CLIENT - For protein component details
# =============================================================================
print("\n" + "="*80)
print("7. TARGET COMPONENT CLIENT")
print("="*80)

components = target_component_client.filter(accession=TEST_UNIPROT_ID)

print(f"\nFound {len(components)} components for {TEST_UNIPROT_ID}")
if components:
    print("\nFirst component:")
    pprint(dict(components[0]))
    
    print("\nAll available fields in component object:")
    pprint(sorted(components[0].keys()))

# =============================================================================
# SUMMARY AND RECOMMENDATIONS
# =============================================================================
print("\n" + "="*80)
print("SUMMARY AND RECOMMENDATIONS")
print("="*80)

print("""
KEY FINDINGS:

1. **ACTIVITY CLIENT**: 
   - Primary source for bioactivity measurements
   - Contains some assay info, but limited
   - Best for: Getting activity values and basic assay identifiers

2. **ASSAY CLIENT**: 
   - ✓ DETAILED assay-level information (organism, cell type, tissue, etc.)
   - ✓ BAO ontology annotations (assay format, label)
   - ✓ Confidence scores
   - ✓ Source database information
   - Should be used to ENRICH activity data with complete assay details

3. **TARGET CLIENT**: 
   - Detailed target information (target type, family, pathways)
   - Protein sequences and classifications
   - Use for: Target characterization and validation

4. **MOLECULE CLIENT**: 
   - Complete molecular properties
   - Structural information
   - Drug development phase
   - Use for: Compound characterization

5. **DOCUMENT CLIENT**: 
   - Full publication metadata
   - Use for: Citation and source tracking

6. **MECHANISM CLIENT**: 
   - Drug mechanism of action
   - Action type (inhibitor, agonist, etc.)
   - Use for: Understanding drug-target interactions

7. **TARGET COMPONENT CLIENT**: 
   - Protein component details
   - Sequence information
   - Use for: Multi-protein target analysis

RECOMMENDED WORKFLOW:
1. Query activities using ACTIVITY CLIENT (as currently done)
2. Extract unique assay_chembl_ids
3. Batch query ASSAY CLIENT to get complete assay details
4. Merge assay details back into activity data
5. Optionally enrich with target/molecule/document details

This approach will give you COMPLETE assay information while maintaining efficiency.
""")

# =============================================================================
# DEMONSTRATION: Proper workflow
# =============================================================================
print("\n" + "="*80)
print("DEMONSTRATION: Recommended Workflow")
print("="*80)

print("\nStep 1: Get activities (current approach)")
activities = activity_client.filter(
    target_organism='Homo sapiens',
    standard_type='IC50',
    target_components__accession=TEST_UNIPROT_ID
)[:10]

print(f"Retrieved {len(activities)} activities")

# Convert to dataframe
activities_df = pd.DataFrame([dict(a) for a in activities])
print(f"\nActivities dataframe shape: {activities_df.shape}")
print(f"Columns: {list(activities_df.columns)}")

print("\nStep 2: Extract unique assay IDs")
unique_assay_ids = activities_df['assay_chembl_id'].unique()
print(f"Found {len(unique_assay_ids)} unique assays")

print("\nStep 3: Batch retrieve detailed assay information")
assay_details_list = []
for assay_id in unique_assay_ids[:3]:  # Demo with first 3
    try:
        assay = assay_client.get(assay_id)
        assay_details_list.append({
            'assay_chembl_id': assay_id,
            'assay_description': assay.get('description'),
            'assay_type': assay.get('assay_type'),
            'assay_type_description': assay.get('assay_type_description'),
            'assay_organism': assay.get('assay_organism'),
            'assay_strain': assay.get('assay_strain'),
            'assay_tissue': assay.get('assay_tissue'),
            'assay_cell_type': assay.get('assay_cell_type'),
            'assay_subcellular_fraction': assay.get('assay_subcellular_fraction'),
            'confidence_score': assay.get('confidence_score'),
            'relationship_type': assay.get('relationship_type'),
            'bao_format': assay.get('bao_format'),
            'bao_label': assay.get('bao_label'),
            'cell_chembl_id': assay.get('cell_chembl_id'),
            'tissue_chembl_id': assay.get('tissue_chembl_id'),
            'src_id': assay.get('src_id'),
            'src_assay_id': assay.get('src_assay_id'),
            'variant_sequence': assay.get('variant_sequence'),
            'target_chembl_id': assay.get('target_chembl_id'),
            'document_chembl_id': assay.get('document_chembl_id')
        })
        print(f"  ✓ Retrieved details for {assay_id}")
    except Exception as e:
        print(f"  ✗ Error for {assay_id}: {e}")

assay_details_df = pd.DataFrame(assay_details_list)
print(f"\nAssay details dataframe shape: {assay_details_df.shape}")

print("\nStep 4: Merge with activity data")
enriched_df = activities_df.merge(
    assay_details_df,
    on='assay_chembl_id',
    how='left',
    suffixes=('_activity', '_assay')
)

print(f"\nEnriched dataframe shape: {enriched_df.shape}")
print(f"Columns: {list(enriched_df.columns)}")

print("\nComparing information richness:")
print(f"  assay_tissue from activity client: {activities_df['assay_tissue'].notna().sum()} non-null values")
print(f"  assay_tissue from assay client: {enriched_df['assay_tissue_assay'].notna().sum()} non-null values")
print(f"  assay_cell_type from activity client: {activities_df['assay_cell_type'].notna().sum()} non-null values")
print(f"  assay_cell_type from assay client: {enriched_df['assay_cell_type_assay'].notna().sum()} non-null values")

# Save sample output
output_file = "data/chembl_API_test/enriched_activities_sample.csv"
enriched_df.to_csv(output_file, index=False)
print(f"\n✓ Sample enriched data saved to {output_file}")

print("\n" + "="*80)
print("Test complete! Check the output for detailed client comparison.")
print("="*80)
