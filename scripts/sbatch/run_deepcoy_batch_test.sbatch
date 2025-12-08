#!/bin/bash
#SBATCH --job-name=deepcoy_batch_test
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_batch_test_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_batch_test_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# =============================================================================
# DeepCoy Batch Test - Process up to 20 UniProt IDs sequentially
#
# This script processes multiple UniProt targets in a single GPU job.
# More efficient than submitting many small jobs for testing.
#
# Usage:
#   sbatch run_deepcoy_batch_test.sbatch
# =============================================================================

# Configuration
PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
DEEPCOY_DIR="${PROJECT_DIR}/external/DeepCoy"
AFFINITY_DIR="${PROJECT_DIR}/data/chembl_affinity"
MAX_TARGETS=20
NUM_DECOYS_PER_ACTIVE=10  # Medium test: 10 decoys per molecule
MAX_MOLECULES=10          # Medium test: 10 molecules per target

# Load conda
module load conda/25.3.0
source activate /scratch1/aoxu/envs/deepcoy

echo "============================================="
echo "DeepCoy Batch Test"
echo "Max targets: ${MAX_TARGETS}"
echo "Molecules per target: ${MAX_MOLECULES}"
echo "Decoys per molecule: ${NUM_DECOYS_PER_ACTIVE}"
echo "============================================="
echo ""

# Find UniProt IDs with filtered parquet files
UNIPROT_IDS=()
for parquet_file in ${AFFINITY_DIR}/uniprot_*/*_filtered.parquet; do
    if [ -f "$parquet_file" ]; then
        dir_name=$(dirname "$parquet_file")
        uid=$(basename "$dir_name" | sed 's/uniprot_//')
        UNIPROT_IDS+=("$uid")
    fi
done

# Limit to MAX_TARGETS
TOTAL_FOUND=${#UNIPROT_IDS[@]}
if [ ${TOTAL_FOUND} -gt ${MAX_TARGETS} ]; then
    UNIPROT_IDS=("${UNIPROT_IDS[@]:0:${MAX_TARGETS}}")
fi

echo "Found ${TOTAL_FOUND} UniProt IDs with filtered parquet"
echo "Processing ${#UNIPROT_IDS[@]} targets"
echo ""

# Track results
SUCCESSFUL=0
FAILED=0
FAILED_IDS=""

# Process each UniProt ID
for i in "${!UNIPROT_IDS[@]}"; do
    UNIPROT_ID="${UNIPROT_IDS[$i]}"
    TARGET_DIR="${AFFINITY_DIR}/uniprot_${UNIPROT_ID}"
    INPUT_PARQUET="${TARGET_DIR}/${UNIPROT_ID}_chembl_activities_filtered.parquet"
    OUTPUT_DIR="${TARGET_DIR}/deepcoy_output_test"
    
    echo "============================================="
    echo "[$((i+1))/${#UNIPROT_IDS[@]}] Processing ${UNIPROT_ID}"
    echo "============================================="
    
    # Check input file exists
    if [ ! -f "$INPUT_PARQUET" ]; then
        echo "ERROR: Input file not found: ${INPUT_PARQUET}"
        FAILED=$((FAILED + 1))
        FAILED_IDS="${FAILED_IDS} ${UNIPROT_ID}"
        continue
    fi
    
    # Create output directory
    mkdir -p ${OUTPUT_DIR}
    
    # Step 1: Convert parquet to DeepCoy JSON format
    echo "Step 1: Converting parquet to JSON..."
    cd ${PROJECT_DIR}
    
    python scripts/prepare_deepcoy_input.py \
        --input_parquet "${INPUT_PARQUET}" \
        --output_json "${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json" \
        --smiles_column canonical_smiles \
        --min_heavy_atoms 10 \
        --max_heavy_atoms 80 \
        --dataset zinc \
        --max_molecules ${MAX_MOLECULES}
    
    if [ ! -f "${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json" ]; then
        echo "ERROR: Failed to create input JSON for ${UNIPROT_ID}"
        FAILED=$((FAILED + 1))
        FAILED_IDS="${FAILED_IDS} ${UNIPROT_ID}"
        continue
    fi
    
    # Check if JSON has any molecules
    NUM_MOLS=$(python -c "import json; print(len(json.load(open('${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json'))))")
    if [ "$NUM_MOLS" -eq 0 ]; then
        echo "WARNING: No valid molecules for ${UNIPROT_ID}, skipping"
        FAILED=$((FAILED + 1))
        FAILED_IDS="${FAILED_IDS} ${UNIPROT_ID}"
        continue
    fi
    echo "Converted ${NUM_MOLS} molecules"
    
    # Step 2: Generate decoys using DeepCoy
    echo "Step 2: Generating decoys..."
    cd ${DEEPCOY_DIR}
    
    python DeepCoy.py \
        --restore models/models/DeepCoy_DUDE_model_e09.pickle \
        --dataset zinc \
        --config "{\"generation\": true, \"number_of_generation_per_valid\": ${NUM_DECOYS_PER_ACTIVE}, \"batch_size\": 1, \"train_file\": \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\", \"valid_file\": \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\", \"output_name\": \"${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt\"}"
    
    if [ -f "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt" ]; then
        NUM_DECOYS=$(wc -l < "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt")
        echo "Generated ${NUM_DECOYS} decoy lines for ${UNIPROT_ID}"
        SUCCESSFUL=$((SUCCESSFUL + 1))
    else
        echo "ERROR: Failed to generate decoys for ${UNIPROT_ID}"
        FAILED=$((FAILED + 1))
        FAILED_IDS="${FAILED_IDS} ${UNIPROT_ID}"
    fi
    
    echo ""
done

# Final summary
echo "============================================="
echo "BATCH TEST COMPLETE"
echo "============================================="
echo "Total targets processed: ${#UNIPROT_IDS[@]}"
echo "Successful: ${SUCCESSFUL}"
echo "Failed: ${FAILED}"
if [ -n "$FAILED_IDS" ]; then
    echo "Failed IDs:${FAILED_IDS}"
fi
echo ""
echo "Output directories: ${AFFINITY_DIR}/uniprot_*/deepcoy_output_test/"
echo "============================================="
