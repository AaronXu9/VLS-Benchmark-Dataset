#!/bin/bash
# =============================================================================
# Dynamic DeepCoy Job Submission - Load Balancing Version (Patched)
#
# Features:
# 1. Scans ALL targets (no fixed ranges)
# 2. Checks for completion (skips generated decoys)
# 3. Randomizes order to minimize collisions between multiple users
# 4. Monitors queue limits
# 5. Atomic locking across multiple users (aoxu + woojinl safe)
# 6. Proper stale-lock handling
#
# Usage:
#   ./submit_deepcoy_dynamic.sh
# =============================================================================

# === CONFIGURATION ===
USER_NAME=$(whoami)
SHARED_AFFINITY_DIR="/project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity"

# User-specific settings
if [ "$USER_NAME" == "aoxu" ]; then
    PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
    SBATCH_SCRIPT="${PROJECT_DIR}/scripts/sbatch/run_deepcoy_single.sbatch"
elif [ "$USER_NAME" == "woojinl" ]; then
    PROJECT_DIR="/scratch1/woojinl/for_Ao/VLS-Benchmark-Dataset"
    SBATCH_SCRIPT="${PROJECT_DIR}/scripts/sbatch/run_deepcoy_single_woojin.sbatch"
else
    echo "Unknown user: $USER_NAME. Please configure paths in script."
    exit 1
fi

MAX_QUEUED=95
CHECK_INTERVAL=60

echo "============================================="
echo "Dynamic DeepCoy Submission"
echo "User:        ${USER_NAME}"
echo "Project Dir: ${PROJECT_DIR}"
echo "Sbatch:      ${SBATCH_SCRIPT}"
echo "============================================="

# Function to count current GPU jobs for THIS user
count_gpu_jobs() {
    squeue -u "$USER_NAME" -p gpu -h 2>/dev/null | wc -l
}

# Function to submit a single job
submit_job() {
    local uid=$1
    result=$(sbatch --export=UNIPROT_ID="${uid}" "${SBATCH_SCRIPT}" 2>&1)
    if echo "$result" | grep -q "Submitted batch job"; then
        echo "[OK] Submitted: ${uid} - $result"
        return 0
    else
        echo "[ERR] Failed:   ${uid} - $result"
        return 1
    fi
}

# Main Loop
while true; do
    # 1. Get list of ALL targets
    # We do this inside the loop to pick up new completions from other users
    echo "Scanning targets..."
    ALL_IDS=()
    PENDING_IDS=()
    COMPLETED_COUNT=0

    # Use glob for ~800 dirs
    for dir in "${SHARED_AFFINITY_DIR}"/uniprot_*; do
        if [ -d "$dir" ]; then
            uid=$(basename "$dir" | sed 's/uniprot_//')

            OUTPUT_DIR="${dir}/deepcoy_output"
            LOCK_FILE="${OUTPUT_DIR}/processing.lock"

            # Completed target: has generated decoys
            if [ -f "${OUTPUT_DIR}/${uid}_generated_decoys.txt" ]; then
                COMPLETED_COUNT=$((COMPLETED_COUNT + 1))

            # Lock present: either active or stale
            elif [ -f "$LOCK_FILE" ]; then
                # Check if lock is stale (> 2 days old)
                if find "$LOCK_FILE" -mtime +2 -print | grep -q .; then
                    echo "[WARN] Found stale lock for ${uid}, clearing and re-queuing..."
                    rm -f "$LOCK_FILE"
                    # Only re-queue if input parquet still exists
                    if [ -f "${dir}/${uid}_chembl_activities_filtered.parquet" ]; then
                        PENDING_IDS+=("$uid")
                    fi
                fi
                # Fresh lock → skip (someone else is working on it)

            else
                # No lock, no completion → valid input?
                if [ -f "${dir}/${uid}_chembl_activities_filtered.parquet" ]; then
                    PENDING_IDS+=("$uid")
                fi
            fi
        fi
    done

    TOTAL_PENDING=${#PENDING_IDS[@]}
    echo "Status: ${COMPLETED_COUNT} completed, ${TOTAL_PENDING} pending (available)."

    if [ "$TOTAL_PENDING" -eq 0 ]; then
        echo "All targets completed! Exiting."
        exit 0
    fi

    # 2. Check Queue (per-user limit)
    CURRENT_JOBS=$(count_gpu_jobs)
    AVAILABLE=$((MAX_QUEUED - CURRENT_JOBS))

    if [ "$AVAILABLE" -le 0 ]; then
        echo "[WAIT] Queue full (${CURRENT_JOBS}/${MAX_QUEUED}). Waiting ${CHECK_INTERVAL}s..."
        sleep "${CHECK_INTERVAL}"
        continue
    fi

    echo "[INFO] Queue has space for ${AVAILABLE} jobs."

    # 3. Shuffle Pending IDs to avoid collisions with other runners
    # We use shuf to randomize the array
    SHUFFLED_IDS=($(shuf -e "${PENDING_IDS[@]}"))

    # 4. Submit Batch
    SUBMITTED_COUNT=0
    for uid in "${SHUFFLED_IDS[@]}"; do
        # Stop if we've filled the available queue slots
        if [ "$SUBMITTED_COUNT" -ge "$AVAILABLE" ]; then
            break
        fi

        LOCK_FILE="${SHARED_AFFINITY_DIR}/uniprot_${uid}/deepcoy_output/processing.lock"
        LOCK_DIR=$(dirname "$LOCK_FILE")
        mkdir -p "$LOCK_DIR"

        # === Atomic lock acquisition ===
        # set -o noclobber; > file  will succeed only if file does not exist.
        # Wrapped in a subshell so noclobber doesn't persist.
        if ( set -o noclobber; > "$LOCK_FILE" ) 2>/dev/null; then
            # We acquired the lock; proceed
            :
        else
            # Another process/user already locked this target; skip
            continue
        fi

        # Now we own the lock for this uid
        if submit_job "$uid"; then
            SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
        else
            # If submission failed, remove lock so others can retry
            rm -f "$LOCK_FILE"
        fi

        # Small delay to be nice to scheduler
        sleep 1
    done

    echo "[INFO] Submitted ${SUBMITTED_COUNT} jobs. Waiting ${CHECK_INTERVAL}s before next cycle..."
    sleep "${CHECK_INTERVAL}"
done