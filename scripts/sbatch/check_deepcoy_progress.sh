#!/bin/bash
# =============================================================================
# Check progress of DeepCoy generation across all targets
#
# Usage:
#   ./check_deepcoy_progress.sh
# =============================================================================

PROJECT_DIR="/scratch1/aoxu/projects/VLS-Benchmark-Dataset"
# Use shared project2 directory accessible by both accounts
AFFINITY_DIR="/project2/katritch_223/aoxu/PLATE-VS/data/chembl_affinity"

echo "============================================="
echo "DeepCoy Generation Progress"
echo "============================================="

# Count total targets
TOTAL=$(ls -d ${AFFINITY_DIR}/uniprot_*/*_filtered.parquet 2>/dev/null | wc -l)

# Count completed (have decoy output)
COMPLETED=$(ls ${AFFINITY_DIR}/uniprot_*/deepcoy_output/*_generated_decoys.txt 2>/dev/null | wc -l)

# Count in progress (have JSON but no decoys)
IN_PROGRESS=0
for json_file in ${AFFINITY_DIR}/uniprot_*/deepcoy_output/molecules_*_actives.json; do
    if [ -f "$json_file" ]; then
        dir=$(dirname "$json_file")
        uid=$(basename $(dirname "$dir") | sed 's/uniprot_//')
        if [ ! -f "${dir}/${uid}_generated_decoys.txt" ]; then
            IN_PROGRESS=$((IN_PROGRESS + 1))
        fi
    fi
done

REMAINING=$((TOTAL - COMPLETED))
PERCENT=$((COMPLETED * 100 / TOTAL))

echo "Total targets:     ${TOTAL}"
echo "Completed:         ${COMPLETED} (${PERCENT}%)"
echo "In progress:       ${IN_PROGRESS}"
echo "Remaining:         ${REMAINING}"
echo ""

# Check running jobs
echo "Running jobs (both accounts):"
squeue -u aoxu -p gpu --format="%.10i %.15j %.8T %.10M" 2>/dev/null | head -20
echo ""
squeue -u woojin -p gpu --format="%.10i %.15j %.8T %.10M" 2>/dev/null | head -20

echo ""
echo "Job counts:"
echo "  aoxu:   $(squeue -u aoxu -p gpu 2>/dev/null | tail -n +2 | wc -l) jobs"
echo "  woojin: $(squeue -u woojin -p gpu 2>/dev/null | tail -n +2 | wc -l) jobs"
echo "============================================="
