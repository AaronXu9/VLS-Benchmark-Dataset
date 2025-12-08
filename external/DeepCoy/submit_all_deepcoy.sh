#!/bin/bash
# =============================================================================
# Submit DeepCoy jobs for all UniProt IDs with filtered parquet files
#
# Usage:
#   ./submit_all_deepcoy.sh           # Submit all jobs
#   ./submit_all_deepcoy.sh --test    # Submit all jobs in test mode
#   ./submit_all_deepcoy.sh --dry-run # Show what would be submitted without submitting
# =============================================================================

PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
AFFINITY_DIR="${PROJECT_DIR}/data/chembl_affinity"
SCRIPT_DIR="${PROJECT_DIR}/external/DeepCoy"

# Parse arguments
TEST_MODE=""
DRY_RUN=false

for arg in "$@"; do
    case $arg in
        --test)
            TEST_MODE="--test"
            ;;
        --dry-run)
            DRY_RUN=true
            ;;
    esac
done

echo "============================================="
echo "DeepCoy Batch Submission"
if [ -n "$TEST_MODE" ]; then
    echo "Mode: TEST"
else
    echo "Mode: FULL"
fi
if [ "$DRY_RUN" = true ]; then
    echo "*** DRY RUN - no jobs will be submitted ***"
fi
echo "============================================="
echo ""

# Find all UniProt directories with filtered parquet files
count=0
submitted=0

for parquet_file in ${AFFINITY_DIR}/uniprot_*/*_filtered.parquet; do
    # Check if file exists (handles case where glob doesn't match)
    if [ ! -f "$parquet_file" ]; then
        continue
    fi
    
    # Extract UniProt ID from path
    dir_name=$(dirname "$parquet_file")
    UNIPROT_ID=$(basename "$dir_name" | sed 's/uniprot_//')
    
    count=$((count + 1))
    
    echo "Found: ${UNIPROT_ID} (${parquet_file})"
    
    if [ "$DRY_RUN" = true ]; then
        echo "  Would submit: sbatch --export=UNIPROT_ID=${UNIPROT_ID} ${SCRIPT_DIR}/run_deepcoy_pipeline.sbatch ${TEST_MODE}"
    else
        # Submit job with UNIPROT_ID as environment variable
        sbatch --export=UNIPROT_ID=${UNIPROT_ID} ${SCRIPT_DIR}/run_deepcoy_pipeline.sbatch ${TEST_MODE}
        submitted=$((submitted + 1))
    fi
    echo ""
done

echo "============================================="
echo "Summary:"
echo "  UniProt IDs found: ${count}"
if [ "$DRY_RUN" = true ]; then
    echo "  Jobs to submit: ${count}"
else
    echo "  Jobs submitted: ${submitted}"
fi
echo "============================================="
