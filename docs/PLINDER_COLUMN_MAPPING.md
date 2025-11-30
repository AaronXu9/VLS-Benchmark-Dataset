# PLINDER Annotation Table Column Mapping

This document maps the filtering criteria to the actual PLINDER annotation table column names based on the official documentation at https://plinder-org.github.io/plinder/

## Quick Reference Table

| Filter Criterion | PLINDER Column Name | Data Type | Description |
|-----------------|-------------------|-----------|-------------|
| **1. Experimental Method** | `entry_determination_method` | str | Experimental method (e.g., "X-RAY DIFFRACTION") |
| **2. Bioassembly** | `system_id` | str | System ID format: `<pdb>__<bioassembly>__<receptor>__<ligand>` |
| **3. Ligand Chain Count** | `system_num_ligand_chains` | int | Number of ligand chains in the system |
| **3. Protein Chain Count** | `system_num_protein_chains` | int | Number of protein chains in the system |
| **4. Covalent Ligands** | `ligand_is_covalent` | list[bool] | Array indicating if each ligand is covalent |
| **5. Ligand Types** | `ligand_is_proper` | list[bool] | Array indicating if each ligand is "proper" (not ion/artifact/cofactor) |
| **6a. Resolution** | `entry_resolution` | float | Resolution in Angstroms |
| **6b. R-factor** | `entry_r` | float | R-factor (R-work) |
| **6c. R-free** | `entry_rfree` | float | R-free factor |
| **6d. Clashscore** | `entry_clashscore` | float | Clashscore from validation |
| **6e. Overall Validation** | `system_pass_validation_criteria` | bool | Pre-computed validation pass/fail |
| **7a. PLIP Interactions** | `system_proper_num_interactions` | int | Number of protein-ligand interactions (proper ligands only) |
| **7b. Pocket Residues** | `system_proper_pocket_num_residues` | int | Number of pocket residues within 6Å (proper ligands only) |
| **8. Molecular Weight** | `system_proper_ligand_max_molecular_weight` | float | Maximum molecular weight of proper ligands (Da) |

## Detailed Column Descriptions

### Entry-Level Columns (PDB-wide)
These columns describe properties of the entire PDB entry:

- **`entry_pdb_id`**: PDB identifier (e.g., "7EEK")
- **`entry_determination_method`**: Experimental method
  - Common values: "X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", "SOLUTION NMR", etc.
  - Filter: Keep only "X-RAY DIFFRACTION" or contains "X-RAY"
  
- **`entry_resolution`**: Resolution in Angstroms (float)
  - Filter: ≤ 3.5 Å
  
- **`entry_rfree`**: R-free factor (float, 0-1 range)
  - Filter: ≤ 0.45
  
- **`entry_r`**: R-factor / R-work (float, 0-1 range)
  - Filter: ≤ 0.4
  - Also check: `entry_rfree - entry_r ≤ 0.05`
  
- **`entry_clashscore`**: Clashscore from validation report
  - Lower is better; typical threshold: < 40

### System-Level Columns
These columns describe a PLINDER system (a specific protein-ligand binding arrangement):

- **`system_id`**: Unique system identifier
  - Format: `<pdb_id>__<bioassembly>__<receptor_chains>__<ligand_chains>`
  - Example: `7eek__1__1.A__1.I` means PDB 7EEK, bioassembly 1, receptor chain 1.A, ligand chain 1.I
  - Extract bioassembly: `system_id.split('__')[1]`
  
- **`system_num_protein_chains`**: Number of protein chains (int)
  - Filter: ≤ 5
  
- **`system_num_ligand_chains`**: Number of ligand chains (int)
  - Filter: ≤ 5

- **`system_pocket_num_residues`**: Total pocket residues (within 6Å of any ligand)
  - Note: This includes ALL ligands (proper + ions + artifacts)

- **`system_num_interactions`**: Total PLIP interactions
  - Note: This includes ALL ligands (proper + ions + artifacts)

### System "Proper" Columns
These columns are calculated using only "proper" ligands (excluding ions, artifacts, cofactors):

- **`system_proper_num_interactions`**: PLIP interactions for proper ligands only
  - Filter: 3-50
  
- **`system_proper_pocket_num_residues`**: Pocket residues for proper ligands only
  - Filter: 5-100
  
- **`system_proper_ligand_max_molecular_weight`**: Max MW of proper ligands (Da)
  - Filter: 200-800

### Ligand-Level Columns (Arrays)
These columns are **arrays** (list[bool], list[float]) because each system can have multiple ligands:

- **`ligand_is_covalent`**: list[bool] - True if ligand is covalently bound
  - Filter: Exclude systems where ANY ligand is covalent
  - Check: `any(ligand_is_covalent)`
  
- **`ligand_is_cofactor`**: list[bool] - True if ligand is a cofactor
  
- **`ligand_is_ion`**: list[bool] - True if ligand is an ion
  
- **`ligand_is_artifact`**: list[bool] - True if ligand is an artifact
  
- **`ligand_is_proper`**: list[bool] - True if ligand is "proper" (not ion/artifact/unwanted type)
  - Filter: Keep systems where at least one ligand is proper
  - Check: `any(ligand_is_proper)`

### Validation Columns (Pre-computed)

- **`system_pass_validation_criteria`**: Boolean - True if system passes PLINDER's crystal quality checks
  - Includes checks for: resolution, R-factors, complete coordinates, no alternate conformations, no clashes, no crystal contacts
  - Recommended to use this in addition to individual metrics
  
- **`system_pass_statistics_criteria`**: Boolean - True if system passes statistics criteria (interaction/pocket counts in acceptable range)

## Implementation Notes

### Handling Array Columns
Many ligand properties are stored as arrays (one value per ligand). Use these patterns:

```python
# Check if ANY ligand meets criterion
has_covalent = df['ligand_is_covalent'].apply(
    lambda x: any(x) if isinstance(x, (list, np.ndarray)) and len(x) > 0 else False
)

# Check if ALL ligands meet criterion
all_proper = df['ligand_is_proper'].apply(
    lambda x: all(x) if isinstance(x, (list, np.ndarray)) and len(x) > 0 else False
)

# Check if at least one ligand meets criterion
has_proper = df['ligand_is_proper'].apply(
    lambda x: any(x) if isinstance(x, (list, np.ndarray)) and len(x) > 0 else False
)
```

### Extracting Bioassembly from system_id
```python
# PLINDER system_id format: pdb__bioassembly__receptor__ligand
bioassembly_nums = df['system_id'].str.split('__').str[1].astype(int)
df_filtered = df[bioassembly_nums == 1]
```

### Using Pre-computed Validation
PLINDER provides `system_pass_validation_criteria` which aggregates many quality checks. It's recommended to use this AND specific filters:

```python
# Apply both specific thresholds AND overall validation
df_filtered = df[
    (df['entry_resolution'] <= 3.5) &
    (df['entry_r'] <= 0.4) &
    (df['entry_rfree'] <= 0.45) &
    (df['system_pass_validation_criteria'] == True)
]
```

## Filter Order and Dependencies

The recommended filtering order (as implemented in the notebook):

1. **X-ray crystallography** - Reduces to ~90% of data
2. **Bioassembly 1** - Reduces to ~60-70% of data
3. **Chain limits** - Removes large complexes
4. **Covalent ligands** - Removes covalent binders
5. **Ligand types** (proper only) - Removes ions/artifacts/cofactors
6. **X-ray validation** - Quality filters (biggest reduction)
7. **Interaction statistics** - Functional relevance
8. **Molecular weight** - Drug-like size range

## Data Source and Version

- **Dataset**: PLINDER v2 (2024-06)
- **Documentation**: https://plinder-org.github.io/plinder/
- **Total Systems**: ~400,000+ protein-ligand systems
- **Total Annotations**: 750+ columns per system
- **File Location**: `data/annotation_table.parquet`

## Known Issues (from PLINDER GitHub)

1. **`entry_release_date`** is INCORRECT - use `query_index` file for correct dates
2. **`ligand_binding_affinity`** columns are DISABLED in current version
3. Some systems may have incorrect nucleic acid receptor classifications

## Example: Complete Filtering Pipeline

```python
import pandas as pd
import numpy as np

# Load data
df = pd.read_parquet('data/annotation_table.parquet')

# Filter 1: X-ray only
df = df[df['entry_determination_method'].str.contains('X-RAY', na=False)]

# Filter 2: Bioassembly 1
bioassembly = df['system_id'].str.split('__').str[1].astype(int)
df = df[bioassembly == 1]

# Filter 3: Chain limits
df = df[(df['system_num_ligand_chains'] <= 5) & 
        (df['system_num_protein_chains'] <= 5)]

# Filter 4: No covalent ligands
has_covalent = df['ligand_is_covalent'].apply(lambda x: any(x) if isinstance(x, (list, np.ndarray)) and len(x) > 0 else False)
df = df[~has_covalent]

# Filter 5: Proper ligands only
has_proper = df['ligand_is_proper'].apply(lambda x: any(x) if isinstance(x, (list, np.ndarray)) and len(x) > 0 else False)
df = df[has_proper]

# Filter 6: X-ray validation
df = df[
    (df['entry_resolution'] <= 3.5) &
    (df['entry_r'] <= 0.4) &
    (df['entry_rfree'] <= 0.45) &
    ((df['entry_rfree'] - df['entry_r']) <= 0.05) &
    (df['system_pass_validation_criteria'] == True)
]

# Filter 7: Interaction statistics (using proper ligand metrics)
df = df[
    (df['system_proper_num_interactions'] >= 3) &
    (df['system_proper_num_interactions'] <= 50) &
    (df['system_proper_pocket_num_residues'] >= 5) &
    (df['system_proper_pocket_num_residues'] <= 100)
]

# Filter 8: Molecular weight
df = df[
    (df['system_proper_ligand_max_molecular_weight'] >= 200) &
    (df['system_proper_ligand_max_molecular_weight'] <= 800)
]

print(f"Final filtered dataset: {len(df):,} systems")
```

## References

1. PLINDER Documentation: https://plinder-org.github.io/plinder/
2. Dataset Tutorial: https://plinder-org.github.io/plinder/tutorial/dataset.html
3. Dataset Reference (Column Descriptions): https://plinder-org.github.io/plinder/dataset.html#annotation-table-target
4. GitHub Repository: https://github.com/plinder-org/plinder
5. Preprint: https://doi.org/10.1101/2024.07.17.603955
