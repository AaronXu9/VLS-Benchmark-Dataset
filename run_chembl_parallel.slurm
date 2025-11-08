#!/bin/bash
#SBATCH --job-name=chembl_parallel
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --output=logs/chembl_retrieval_%j.out
#SBATCH --error=logs/chembl_retrieval_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<YOUR_EMAIL>

# USC CARC SLURM Job Script for Parallel ChEMBL Retrieval
# 
# Usage:
#   sbatch run_chembl_parallel.slurm
#
# To customize resources:
#   sbatch --cpus-per-task=32 --mem=64GB run_chembl_parallel.slurm

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Number of CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE MB"
echo "Start Time: $(date)"
echo "=================================================="

# Load required modules
module purge
module load gcc/11.3.0
module load python/3.11.3

# Activate your Python environment
# Option 1: If using conda
source ~/.bashrc
conda activate plot  # or your environment name

# Option 2: If using venv
# source /path/to/your/venv/bin/activate

# Set working directory
cd $SLURM_SUBMIT_DIR
echo "Working directory: $(pwd)"

# Create logs directory if it doesn't exist
mkdir -p logs

# Print Python and package info
echo ""
echo "Python version:"
python --version
echo ""
echo "Installed packages:"
pip list | grep -E "pandas|numpy|chembl|requests"
echo ""

# Run the parallel retrieval script
echo "=================================================="
echo "Starting ChEMBL parallel retrieval..."
echo "=================================================="

python retrieve_chembl_parallel.py \
    --input data/annotation_table_with_uniprot.parquet \
    --output data/chembl_parallel \
    --final-output data/chembl_activities_merged.parquet \
    --workers $SLURM_CPUS_PER_TASK \
    --max-activities 1000 \
    --standard-types IC50 Ki Kd EC50

EXIT_CODE=$?

echo ""
echo "=================================================="
echo "Job finished with exit code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=================================================="

# Print summary of output files
if [ -f data/chembl_parallel/retrieval_summary.json ]; then
    echo ""
    echo "Retrieval Summary:"
    python -c "import json; data=json.load(open('data/chembl_parallel/retrieval_summary.json')); print(f\"  Total UniProt IDs: {data['total_uniprot_ids']}\"); print(f\"  Successful: {data['successful']}\"); print(f\"  Failed: {data['failed']}\")"
fi

exit $EXIT_CODE
