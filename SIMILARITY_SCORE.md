

# Usage 
```python
from pli_similarity_scorer import PLISimilarityScorer, load_annotation_table
from pathlib import Path

# Load data
annotation_table = load_annotation_table(Path('data/matched_annotation_table.parquet'))

# Create scorer
scorer = PLISimilarityScorer(
    systems_dir=Path('/mnt/katritch_lab2/aoxu/2024-06/v2/systems'),
    annotation_table=annotation_table,
)

# Score a pair
score = scorer.score_pair('4grb__1__1.A__1.C', '6gra__1__1.A__1.B')
print(f"Tanimoto: {score.tanimoto:.4f}")
print(f"SuCOS (protein-aligned): {score.sucos_protein_aligned:.4f}")
print(f"SuCOS (ligand-aligned): {score.sucos_ligand_aligned:.4f}")

# Batch scoring
pairs = [('sys1', 'sys2'), ('sys3', 'sys4')]
df = scorer.score_batch(pairs)
