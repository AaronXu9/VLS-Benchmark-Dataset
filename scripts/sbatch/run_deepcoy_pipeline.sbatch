#!/bin/bash
#SBATCH --job-name=deepcoy
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_%j.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1

# =============================================================================
# DeepCoy Pipeline for UniProt targets
# 
# This script:
# 1. Converts ChEMBL affinity parquet data to DeepCoy JSON input format
# 2. Generates decoy molecules using pretrained DeepCoy model
# 3. Evaluates generated decoys for quality
#
# Usage:
#   Single target:  sbatch --export=UNIPROT_ID=P10632 run_deepcoy_pipeline.sbatch
#   Test mode:      sbatch --export=UNIPROT_ID=P10632 run_deepcoy_pipeline.sbatch --test
#   Batch:          ./submit_all_deepcoy.sh
# =============================================================================

# Parse command line arguments
TEST_MODE=false
for arg in "$@"; do
    case $arg in
        --test)
            TEST_MODE=true
            shift
            ;;
    esac
done

# Configuration - UNIPROT_ID can be set via --export or defaults to P10632
UNIPROT_ID="${UNIPROT_ID:-P10632}"
PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
DEEPCOY_DIR="${PROJECT_DIR}/external/DeepCoy"
AFFINITY_DIR="${PROJECT_DIR}/data/chembl_affinity/uniprot_${UNIPROT_ID}"
INPUT_PARQUET="${AFFINITY_DIR}/${UNIPROT_ID}_chembl_activities_filtered.parquet"

# Test mode settings: fewer molecules, fewer decoys, skip evaluation
if [ "$TEST_MODE" = true ]; then
    echo ">>> RUNNING IN TEST MODE <<<"
    OUTPUT_DIR="${AFFINITY_DIR}/deepcoy_output_test"
    NUM_DECOYS_PER_ACTIVE=5
    MAX_MOLECULES=3
    SKIP_EVALUATION=true
else
    OUTPUT_DIR="${AFFINITY_DIR}/deepcoy_output"
    NUM_DECOYS_PER_ACTIVE=50
    MAX_MOLECULES=0  # 0 means no limit
    SKIP_EVALUATION=false
fi

# Load conda
module load conda/25.3.0

# Activate the deepcoy environment
source activate /scratch1/aoxu/envs/deepcoy

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "============================================="
echo "DeepCoy Pipeline for UniProt ${UNIPROT_ID}"
if [ "$TEST_MODE" = true ]; then
    echo ">>> TEST MODE: ${MAX_MOLECULES} molecules, ${NUM_DECOYS_PER_ACTIVE} decoys each <<<"
fi
echo "============================================="
echo "Input: ${INPUT_PARQUET}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# =============================================================================
# Step 1: Convert parquet to DeepCoy JSON format
# =============================================================================
echo "Step 1: Converting parquet to DeepCoy JSON format..."

cd ${PROJECT_DIR}

# Build command with optional --max_molecules for test mode
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
    echo "ERROR: Failed to create input JSON file"
    exit 1
fi

echo "Created: ${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json"
echo ""

# =============================================================================
# Step 2: Generate decoys using DeepCoy
# =============================================================================
echo "Step 2: Generating decoys with DeepCoy..."

cd ${DEEPCOY_DIR}

python DeepCoy.py \
    --restore models/models/DeepCoy_DUDE_model_e09.pickle \
    --dataset zinc \
    --config "{\"generation\": true, \"number_of_generation_per_valid\": ${NUM_DECOYS_PER_ACTIVE}, \"batch_size\": 1, \"train_file\": \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\", \"valid_file\": \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\", \"output_name\": \"${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt\"}"

if [ ! -f "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt" ]; then
    echo "ERROR: Failed to generate decoys"
    exit 1
fi

echo "Generated decoys: ${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt"
echo "Number of generated lines:"
wc -l "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt"
echo ""

# =============================================================================
# Step 3: Evaluate generated decoys (skip in test mode)
# =============================================================================
if [ "$SKIP_EVALUATION" = true ]; then
    echo "Step 3: SKIPPED (test mode)"
else
    echo "Step 3: Evaluating generated decoys..."

    mkdir -p ${OUTPUT_DIR}/evaluation

    python evaluation/select_and_evaluate_decoys.py \
        --data_path "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt" \
        --output_path "${OUTPUT_DIR}/evaluation/" \
        --dataset_name dude \
        --num_decoys_per_active 50 \
        --min_num_candidates 100 \
        --num_cores 4
fi

echo ""
echo "============================================="
echo "Pipeline Complete!"
echo "============================================="
echo ""
echo "Output files (in ${AFFINITY_DIR}/deepcoy_output/):"
echo "  - Input molecules: ${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json"
echo "  - Generated decoys: ${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt"
echo "  - Evaluation results: ${OUTPUT_DIR}/evaluation/results.csv"
echo "  - Selected decoys: ${OUTPUT_DIR}/evaluation/${UNIPROT_ID}_generated_decoys-selected.smi"
echo ""

if [ -f "${OUTPUT_DIR}/evaluation/results.csv" ]; then
    echo "Evaluation Results:"
    cat "${OUTPUT_DIR}/evaluation/results.csv"
fi
