#!/bin/bash
#SBATCH --job-name=deepcoy_gen
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/logs/deepcoy_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/logs/deepcoy_%j.err

# Load conda module
module purge
module load conda/25.3.0

# Activate DeepCoy environment
source activate /scratch1/aoxu/envs/deepcoy

# Set working directory
cd /scratch1/aoxu/projects/VLS-Benchmark-Dataset/external/DeepCoy

# Create output directory if it doesn't exist
mkdir -p output

# Run DeepCoy generation with pretrained DUDE model
# Using the validation set as input molecules to generate decoys
python DeepCoy.py \
    --restore models/models/DeepCoy_DUDE_model_e09.pickle \
    --dataset zinc \
    --config '{
        "generation": true,
        "number_of_generation_per_valid": 100,
        "batch_size": 1,
        "train_file": "data/molecules_zinc_dekois_valid.json",
        "valid_file": "data/molecules_zinc_dekois_valid.json",
        "output_name": "output/DeepCoy_generated_decoys_test.txt"
    }'

echo "DeepCoy generation completed!"
