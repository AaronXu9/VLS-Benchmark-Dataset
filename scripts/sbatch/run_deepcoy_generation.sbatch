#!/bin/bash
#SBATCH --job-name=deepcoy_setup
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_setup_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_setup_%j.err
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# Load conda
module load conda/25.3.0

# Activate the deepcoy environment
source activate /scratch1/aoxu/envs/deepcoy

# Install missing dependencies
pip install docopt

# Change to data directory and prepare data
cd /scratch1/aoxu/projects/VLS-Benchmark-Dataset/external/DeepCoy/data
python prepare_data.py

echo "Data preparation completed!"
echo "Listing generated files:"
ls -la *.json 2>/dev/null || echo "No JSON files generated yet"