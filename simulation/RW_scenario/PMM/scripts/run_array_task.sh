#!/usr/bin/env bash
set -euo pipefail

pmm_dir=$(cd "$(dirname "$TASK_FILE")/.." && pwd)
script_dir="$pmm_dir/scripts"
module load R/4.5
export R_LIBS_USER="${R_LIBS_USER:-$HOME/MI/Rlib}"
unset TMPDIR TMP TEMP
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
echo "$(timestamp) | job_id=${SLURM_JOB_ID:-local} task_id=${SLURM_ARRAY_TASK_ID:-local} state=START"
runtime_dir=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/rw_pmm_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-local}_XXXXXX")
trap 'rm -rf "$runtime_dir"' EXIT
cd "$runtime_dir"
Rscript "$script_dir/simulate_batch.R" \
  "$TASK_FILE" "$SLURM_ARRAY_TASK_ID" "$RAW_DIR"
echo "$(timestamp) | job_id=${SLURM_JOB_ID:-local} task_id=${SLURM_ARRAY_TASK_ID:-local} state=COMPLETE"
