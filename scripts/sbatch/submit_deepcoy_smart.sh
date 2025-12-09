#!/bin/bash
# =============================================================================
# Smart DeepCoy Job Submission with Queue Management - AOXU ACCOUNT
#
# This script monitors queue limits and submits jobs only when there's room.
# It will keep running until all targets are submitted.
#
# Usage:
#   ./submit_deepcoy_smart.sh                    # Submit for current user
#   ./submit_deepcoy_smart.sh --start 0 --end 418    # Submit specific range
#   ./submit_deepcoy_smart.sh --dry-run          # Preview without submitting
#   nohup ./submit_deepcoy_smart.sh > submit.log 2>&1 &  # Run in background
# =============================================================================

# === AOXU ACCOUNT PATHS ===
PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
# Use shared project2 directory accessible by both accounts
AFFINITY_DIR="/project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity"
SBATCH_SCRIPT="${PROJECT_DIR}/scripts/sbatch/run_deepcoy_single.sbatch"

# Queue limits
MAX_QUEUED=95       # Stay under 100 limit with buffer
CHECK_INTERVAL=60   # Seconds between queue checks

# Parse arguments
DRY_RUN=false
START_IDX=-1
END_IDX=-1

while [[ $# -gt 0 ]]; do
    case $1 in
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --start)
            START_IDX=$2
            shift 2
            ;;
        --end)
            END_IDX=$2
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

# Get all UniProt IDs with filtered parquet
ALL_IDS=()
for parquet_file in ${AFFINITY_DIR}/uniprot_*/*_filtered.parquet; do
    if [ -f "$parquet_file" ]; then
        dir_name=$(dirname "$parquet_file")
        uid=$(basename "$dir_name" | sed 's/uniprot_//')
        ALL_IDS+=("$uid")
    fi
done

# Sort for consistent ordering
IFS=$'\n' SORTED_IDS=($(sort <<<"${ALL_IDS[*]}")); unset IFS
TOTAL=${#SORTED_IDS[@]}

# Set default range for aoxu (first half: 0 to TOTAL/2-1)
if [ $START_IDX -eq -1 ]; then
    START_IDX=0
fi
if [ $END_IDX -eq -1 ]; then
    END_IDX=$((TOTAL / 2 - 1))
fi

echo "============================================="
echo "Smart DeepCoy Job Submission - AOXU"
echo "============================================="
echo "Total targets available: ${TOTAL}"
echo "Processing range: ${START_IDX} to ${END_IDX}"
echo "Max queued jobs: ${MAX_QUEUED}"
echo "Check interval: ${CHECK_INTERVAL}s"
if [ "$DRY_RUN" = true ]; then
    echo "*** DRY RUN MODE ***"
fi
echo "============================================="
echo ""

# Function to count current GPU jobs
count_gpu_jobs() {
    squeue -u $USER -p gpu -h 2>/dev/null | wc -l
}

# Function to submit a single job
submit_job() {
    local uid=$1
    if [ "$DRY_RUN" = true ]; then
        echo "[DRY] Would submit: ${uid}"
        return 0
    else
        result=$(sbatch --export=UNIPROT_ID=${uid} ${SBATCH_SCRIPT} 2>&1)
        if echo "$result" | grep -q "Submitted batch job"; then
            echo "[OK] Submitted: ${uid} - $result"
            return 0
        else
            echo "[ERR] Failed: ${uid} - $result"
            return 1
        fi
    fi
}

# Build list of targets to submit (skip completed ones)
PENDING_IDS=()
SKIPPED=0
for i in $(seq ${START_IDX} ${END_IDX}); do
    UNIPROT_ID="${SORTED_IDS[$i]}"
    OUTPUT_DIR="${AFFINITY_DIR}/uniprot_${UNIPROT_ID}/deepcoy_output"
    
    # Skip if decoys already generated
    if [ -f "${OUTPUT_DIR}/${UNIPROT_ID}_generated_decoys.txt" ]; then
        SKIPPED=$((SKIPPED + 1))
        continue  # Skip completed
    fi
    PENDING_IDS+=("$UNIPROT_ID")
done

echo "Targets to submit: ${#PENDING_IDS[@]}"
echo "Already completed (skipped): ${SKIPPED}"
echo ""

# Submit jobs with queue management
SUBMITTED=0
FAILED=0
IDX=0

while [ $IDX -lt ${#PENDING_IDS[@]} ]; do
    # Check current queue size
    CURRENT_JOBS=$(count_gpu_jobs)
    AVAILABLE=$((MAX_QUEUED - CURRENT_JOBS))
    
    if [ $AVAILABLE -le 0 ]; then
        echo "[WAIT] Queue full (${CURRENT_JOBS}/${MAX_QUEUED}). Waiting ${CHECK_INTERVAL}s..."
        sleep ${CHECK_INTERVAL}
        continue
    fi
    
    # Submit up to AVAILABLE jobs
    BATCH_COUNT=0
    while [ $IDX -lt ${#PENDING_IDS[@]} ] && [ $BATCH_COUNT -lt $AVAILABLE ]; do
        UNIPROT_ID="${PENDING_IDS[$IDX]}"
        
        if submit_job "$UNIPROT_ID"; then
            SUBMITTED=$((SUBMITTED + 1))
        else
            FAILED=$((FAILED + 1))
        fi
        
        IDX=$((IDX + 1))
        BATCH_COUNT=$((BATCH_COUNT + 1))
        
        # Small delay between submissions
        sleep 0.5
    done
    
    echo "[INFO] Progress: ${SUBMITTED}/${#PENDING_IDS[@]} submitted, ${FAILED} failed"
    
    # If more jobs to submit, wait before next batch
    if [ $IDX -lt ${#PENDING_IDS[@]} ]; then
        echo "[WAIT] Waiting ${CHECK_INTERVAL}s before next batch..."
        sleep ${CHECK_INTERVAL}
    fi
done

echo ""
echo "============================================="
echo "SUBMISSION COMPLETE - AOXU"
echo "============================================="
echo "Total submitted: ${SUBMITTED}"
echo "Failed: ${FAILED}"
echo "Skipped (already done): ${SKIPPED}"
echo "============================================="
