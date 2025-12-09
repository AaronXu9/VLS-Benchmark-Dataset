#!/bin/bash
# =============================================================================
# Submit DeepCoy jobs for account: woojin (second half of targets)
#
# Usage:
#   ./submit_deepcoy_woojin.sh           # Submit jobs
#   ./submit_deepcoy_woojin.sh --dry-run # Preview without submitting
# =============================================================================

PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
# Use shared project2 directory accessible by both accounts
AFFINITY_DIR="/project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity"
SBATCH_SCRIPT="${PROJECT_DIR}/scripts/sbatch/run_deepcoy_single_woojin.sbatch"

# Account settings
ACCOUNT="woojin"
MAX_QUEUED=100  # GPU queue limit

# Parse arguments
DRY_RUN=false
if [ "$1" == "--dry-run" ]; then
    DRY_RUN=true
fi

# Get all UniProt IDs with filtered parquet
ALL_IDS=()
for parquet_file in ${AFFINITY_DIR}/uniprot_*/*_filtered.parquet; do
    if [ -f "$parquet_file" ]; then
        dir_name=$(dirname "$parquet_file")
        uid=$(basename "$dir_name" | sed 's/uniprot_//')
        ALL_IDS+=("$uid")
    fi
done

# Sort for consistent splitting
IFS=$'\n' SORTED_IDS=($(sort <<<"${ALL_IDS[*]}")); unset IFS

TOTAL=${#SORTED_IDS[@]}
HALF=$((TOTAL / 2))

# Account woojin gets second half (indices HALF to TOTAL-1)
START_IDX=${HALF}
END_IDX=$((TOTAL - 1))
COUNT=$((TOTAL - HALF))

echo "============================================="
echo "DeepCoy Job Submission - Account: ${ACCOUNT}"
echo "============================================="
echo "Total targets: ${TOTAL}"
echo "This account processes: ${START_IDX} to ${END_IDX} (${COUNT} targets)"
echo "Max queued jobs: ${MAX_QUEUED}"
if [ "$DRY_RUN" = true ]; then
    echo "*** DRY RUN MODE ***"
fi
echo ""

# Submit jobs
SUBMITTED=0
SKIPPED=0

for i in $(seq ${START_IDX} ${END_IDX}); do
    UNIPROT_ID="${SORTED_IDS[$i]}"
    OUTPUT_DIR="${AFFINITY_DIR}/uniprot_${UNIPROT_ID}/deepcoy_output"
    
    # Skip if already completed
    if [ -f "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt" ]; then
        echo "SKIP (done): ${UNIPROT_ID}"
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    
    if [ "$DRY_RUN" = true ]; then
        echo "Would submit: ${UNIPROT_ID}"
    else
        sbatch --export=UNIPROT_ID=${UNIPROT_ID} ${SBATCH_SCRIPT}
        echo "Submitted: ${UNIPROT_ID}"
    fi
    
    SUBMITTED=$((SUBMITTED + 1))
    
    # Pause every 50 jobs to avoid overwhelming scheduler
    if [ $((SUBMITTED % 50)) -eq 0 ] && [ "$DRY_RUN" = false ]; then
        echo "Pausing 5s after ${SUBMITTED} submissions..."
        sleep 5
    fi
done

echo ""
echo "============================================="
echo "Summary for ${ACCOUNT}:"
echo "  Submitted: ${SUBMITTED}"
echo "  Skipped (already done): ${SKIPPED}"
echo "============================================="
