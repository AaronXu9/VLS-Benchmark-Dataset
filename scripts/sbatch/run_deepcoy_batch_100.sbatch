#!/bin/bash
#SBATCH --job-name=deepcoy_batch_100
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_batch_100_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_batch_100_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-task=p100:1
#SBATCH --mem=64G
#SBATCH --partition=gpu

# =============================================================================
# DeepCoy Batch - Generate 100 decoys per active for 10 targets
#
# This is a production-quality test with enough decoys for meaningful evaluation.
# Settings:
#   - 10 UniProt targets
#   - All actives per target (no limit)
#   - 100 decoys per active molecule
#
# Usage:
#   sbatch run_deepcoy_batch_100.sbatch
# =============================================================================

# Configuration
PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
DEEPCOY_DIR="${PROJECT_DIR}/external/DeepCoy"
AFFINITY_DIR="${PROJECT_DIR}/data/chembl_affinity"
MAX_TARGETS=10
NUM_DECOYS_PER_ACTIVE=100  # 100 decoys per molecule for good selection
MAX_MOLECULES=0             # 0 = no limit, use all actives

# Output directory suffix
OUTPUT_SUFFIX="deepcoy_output_100"

# Load conda
module load conda/25.3.0
source activate /scratch1/aoxu/envs/deepcoy

echo "============================================="
echo "DeepCoy Batch - 100 Decoys per Active"
echo "Max targets: ${MAX_TARGETS}"
echo "Molecules per target: ALL (no limit)"
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
    OUTPUT_DIR="${TARGET_DIR}/${OUTPUT_SUFFIX}"
    
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
    
    # Build command - no --max_molecules means use all
    CONVERT_CMD="python scripts/prepare_deepcoy_input.py \
        --input_parquet \"${INPUT_PARQUET}\" \
        --output_json \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\" \
        --smiles_column canonical_smiles \
        --min_heavy_atoms 10 \
        --max_heavy_atoms 80 \
        --dataset zinc"
    
    if [ "$MAX_MOLECULES" -gt 0 ]; then
        CONVERT_CMD="${CONVERT_CMD} --max_molecules ${MAX_MOLECULES}"
    fi
    
    eval ${CONVERT_CMD}
    
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
    echo "Step 2: Generating ${NUM_DECOYS_PER_ACTIVE} decoys per active..."
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
echo "BATCH GENERATION COMPLETE"
echo "============================================="
echo "Total targets processed: ${#UNIPROT_IDS[@]}"
echo "Successful: ${SUCCESSFUL}"
echo "Failed: ${FAILED}"
if [ -n "$FAILED_IDS" ]; then
    echo "Failed IDs:${FAILED_IDS}"
fi
echo ""
echo "Output directories: ${AFFINITY_DIR}/uniprot_*/${OUTPUT_SUFFIX}/"
echo ""
echo "Next step: Run evaluation with:"
echo "  sbatch run_deepcoy_batch_100_eval.sbatch"
echo "============================================="
