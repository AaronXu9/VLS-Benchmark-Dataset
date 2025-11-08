#!/bin/bash
#SBATCH --job-name=uniprot_extract
#SBATCH --account=katritch_223
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16GB
#SBATCH --time=06:00:00
#SBATCH --output=logs/uniprot_extract_%j.out
#SBATCH --error=logs/uniprot_extract_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aoxu@usc.edu

# Load required modules
module purge
module load python/3.11.3
module load gcc/11.3.0

# Activate conda environment
source ~/.bashrc
conda activate PoseBench 

# Create output directories if they don't exist
mkdir -p logs
mkdir -p data

# Print job information
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Start Time: $(date)"
echo "=========================================="

# Run the UniProt extraction script
python extract_uniprot_parallel.py \
    --input data/annotation_table_filtered.parquet \
    --output data/annotation_table_with_uniprot.parquet \
    --workers 16 \
    --batch-size 20

# Check exit status
EXIT_CODE=$?
echo "=========================================="
echo "End Time: $(date)"
echo "Exit Code: $EXIT_CODE"

if [ $EXIT_CODE -eq 0 ]; then
    echo "UniProt extraction completed successfully!"
    
    # Print summary statistics
    echo ""
    echo "Output file statistics:"
    python -c "
import pandas as pd
df = pd.read_parquet('data/annotation_table_with_uniprot.parquet')
print(f'Total systems: {len(df)}')
print(f'Systems with UniProt IDs: {df[\"uniprot_ids_str\"].astype(bool).sum()}')
print(f'Unique UniProt IDs: {df[\"uniprot_ids_str\"].str.split(\",\").explode().nunique()}')
"
else
    echo "UniProt extraction failed with exit code $EXIT_CODE"
    exit $EXIT_CODE
fi
