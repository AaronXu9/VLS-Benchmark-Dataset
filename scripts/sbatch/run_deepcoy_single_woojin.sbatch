#!/bin/bash
#SBATCH --job-name=dcoy_%j
#SBATCH --output=/scratch1/woojinl/for_Ao/VLS-Benchmark-Dataset/slurm/deepcoy_%j.out
#SBATCH --error=/scratch1/woojinl/for_Ao/VLS-Benchmark-Dataset/slurm/deepcoy_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=gpu
#SBATCH --gpus-per-task=p100:1

# =============================================================================
# DeepCoy Single Target Job - WOOJIN ACCOUNT
#
# Time estimation: ~0.573s per decoy
#   - 20 decoys x 100 actives = 2000 decoys = ~19 min
#   - 20 decoys x 500 actives = 10000 decoys = ~96 min
#   - 4 hours should cover most targets with buffer
# 
# Usage:
#   sbatch --export=UNIPROT_ID=P12345 run_deepcoy_single_woojin.sbatch
# =============================================================================

# Check UNIPROT_ID is set
if [ -z "$UNIPROT_ID" ]; then
    echo "ERROR: UNIPROT_ID not set. Use --export=UNIPROT_ID=..."
    exit 1
fi

# === WOOJIN ACCOUNT PATHS ===
PROJECT_DIR="/scratch1/woojinl/for_Ao/VLS-Benchmark-Dataset"
DEEPCOY_DIR="${PROJECT_DIR}/external/DeepCoy"
# Use shared project2 directory accessible by both accounts
AFFINITY_DIR="/project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity"
NUM_DECOYS_PER_ACTIVE=5
OUTPUT_SUFFIX="deepcoy_output"

TARGET_DIR="${AFFINITY_DIR}/uniprot_${UNIPROT_ID}"
INPUT_PARQUET="${TARGET_DIR}/${UNIPROT_ID}_chembl_activities_filtered.parquet"
OUTPUT_DIR="${TARGET_DIR}/${OUTPUT_SUFFIX}"

# Load environment - WOOJIN's conda
module load conda/25.3.0
source activate /scratch1/woojinl/envs/deepcoy

# Fix protobuf compatibility with TensorFlow 1.x
export PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

echo "============================================="
echo "DeepCoy Generation: ${UNIPROT_ID}"
echo "Account: woojin"
echo "Decoys per active: ${NUM_DECOYS_PER_ACTIVE}"
echo "============================================="

# Check input exists
if [ ! -f "$INPUT_PARQUET" ]; then
    echo "ERROR: Input not found: ${INPUT_PARQUET}"
    exit 1
fi

mkdir -p ${OUTPUT_DIR}

# Step 1: Convert parquet to JSON
echo "Step 1: Converting parquet to JSON..."
cd ${PROJECT_DIR}

python scripts/prepare_deepcoy_input.py \
    --input_parquet "${INPUT_PARQUET}" \
    --output_json "${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json" \
    --smiles_column canonical_smiles \
    --min_heavy_atoms 10 \
    --max_heavy_atoms 80 \
    --dataset zinc

if [ ! -f "${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json" ]; then
    echo "ERROR: Failed to create JSON"
    exit 1
fi

NUM_MOLS=$(python -c "import json; print(len(json.load(open('${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json'))))")
echo "Molecules: ${NUM_MOLS}"

if [ "$NUM_MOLS" -eq 0 ]; then
    echo "WARNING: No valid molecules, skipping"
    exit 0
fi

# Step 2: Generate decoys
echo "Step 2: Generating decoys..."
cd ${DEEPCOY_DIR}

python DeepCoy.py \
    --restore models/models/DeepCoy_DUDE_model_e09.pickle \
    --dataset zinc \
    --config "{\"generation\": true, \"number_of_generation_per_valid\": ${NUM_DECOYS_PER_ACTIVE}, \"batch_size\": 1, \"train_file\": \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\", \"valid_file\": \"${OUTPUT_DIR}/molecules_${UNIPROT_ID}_actives.json\", \"output_name\": \"${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt\"}"

if [ -f "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt" ]; then
    NUM_DECOYS=$(wc -l < "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt")
    echo "SUCCESS: Generated ${NUM_DECOYS} decoys for ${UNIPROT_ID}"
    # Remove lock file on success
    rm -f "${OUTPUT_DIR}/processing.lock"
else
    echo "ERROR: Failed to generate decoys"
    # Remove lock file on failure too (so it can be retried)
    rm -f "${OUTPUT_DIR}/processing.lock"
    exit 1
fi
