#!/bin/bash
#SBATCH --job-name=deepcoy_benchmark
#SBATCH --output=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_benchmark_%j.out
#SBATCH --error=/scratch1/aoxu/projects/VLS-Benchmark-Dataset/slurm/deepcoy_benchmark_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gpus-per-task=p100:1
#SBATCH --mem=16G
#SBATCH --partition=gpu

# =============================================================================
# DeepCoy Benchmark - Measure time per decoy generation
# =============================================================================

PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
DEEPCOY_DIR="${PROJECT_DIR}/external/DeepCoy"

module load conda/25.3.0
source activate /scratch1/aoxu/envs/deepcoy

cd ${DEEPCOY_DIR}

echo "============================================="
echo "DeepCoy Benchmark"
echo "============================================="

# Use existing test input (3 molecules)
TEST_INPUT="${PROJECT_DIR}/data/chembl_affinity/uniprot_P10632/deepcoy_output_test/molecules_P10632_actives.json"

if [ ! -f "$TEST_INPUT" ]; then
    echo "Creating test input..."
    cd ${PROJECT_DIR}
    python scripts/prepare_deepcoy_input.py \
        --input_parquet data/chembl_affinity/uniprot_P10632/P10632_chembl_activities_filtered.parquet \
        --output_json ${TEST_INPUT} \
        --max_molecules 5
    cd ${DEEPCOY_DIR}
fi

NUM_MOLS=$(python -c "import json; print(len(json.load(open('${TEST_INPUT}'))))")
echo "Test molecules: ${NUM_MOLS}"

# Benchmark different decoy counts
for NUM_DECOYS in 10 50 100; do
    echo ""
    echo "--- Generating ${NUM_DECOYS} decoys per molecule ---"
    
    START_TIME=$(date +%s.%N)
    
    python DeepCoy.py \
        --restore models/models/DeepCoy_DUDE_model_e09.pickle \
        --dataset zinc \
        --config "{\"generation\": true, \"number_of_generation_per_valid\": ${NUM_DECOYS}, \"batch_size\": 1, \"train_file\": \"${TEST_INPUT}\", \"valid_file\": \"${TEST_INPUT}\", \"output_name\": \"benchmark_output.txt\"}" \
        2>&1 | tail -5
    
    END_TIME=$(date +%s.%N)
    ELAPSED=$(echo "$END_TIME - $START_TIME" | bc)
    
    TOTAL_DECOYS=$((NUM_MOLS * NUM_DECOYS))
    TIME_PER_DECOY=$(echo "scale=4; $ELAPSED / $TOTAL_DECOYS" | bc)
    TIME_PER_MOL=$(echo "scale=4; $ELAPSED / $NUM_MOLS" | bc)
    
    echo "Total time: ${ELAPSED}s"
    echo "Time per molecule: ${TIME_PER_MOL}s"
    echo "Time per decoy: ${TIME_PER_DECOY}s"
    echo "Total decoys generated: ${TOTAL_DECOYS}"
done

# Clean up
rm -f benchmark_output.txt

echo ""
echo "============================================="
echo "ESTIMATION CALCULATOR"
echo "============================================="
echo ""
echo "Use these formulas to estimate your job time:"
echo "  Total time = (num_actives) x (num_decoys) x (time_per_decoy)"
echo ""
echo "Example for 100 actives x 100 decoys:"
echo "  If time_per_decoy = 0.5s: 100 x 100 x 0.5 = 5000s = 1.4 hours"
echo "  If time_per_decoy = 1.0s: 100 x 100 x 1.0 = 10000s = 2.8 hours"
echo "============================================="