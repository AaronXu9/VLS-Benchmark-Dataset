#!/bin/bash
#SBATCH --job-name=deepcoy_eval
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_eval_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_eval_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=main

# Load conda
module load conda/25.3.0

# Activate the deepcoy environment
source activate /scratch1/aoxu/envs/deepcoy

# Change to DeepCoy directory
cd /scratch1/aoxu/projects/VLS-Benchmark-Dataset/external/DeepCoy

# Create output directory for evaluation results
mkdir -p evaluation_output

echo "Starting DeepCoy evaluation..."
echo "Input file: output/P38-alpha_generated_decoys.txt"
echo "Output directory: evaluation_output/"

# Run evaluation script
python evaluation/select_and_evaluate_decoys.py \
    --data_path output/P38-alpha_generated_decoys.txt \
    --output_path evaluation_output/ \
    --dataset_name dude \
    --num_decoys_per_active 50 \
    --min_num_candidates 100 \
    --num_cores 4

echo ""
echo "Evaluation completed!"
echo ""
echo "Results summary:"
cat evaluation_output/results.csv
echo ""
echo "Selected decoys saved to: evaluation_output/P38-alpha_generated_decoys-selected.smi"
wc -l evaluation_output/P38-alpha_generated_decoys-selected.smi
