#!/usr/bin/env bash
set -euo pipefail

pmm_dir=$(cd "$(dirname "$TASK_FILE")" && pwd)
script_dir="$pmm_dir/scripts"
if [[ -n "${RSCRIPT_BIN:-}" ]]; then
  rscript="$RSCRIPT_BIN"
elif command -v Rscript >/dev/null 2>&1; then
  rscript="Rscript"
else
  module load R/4.5
  rscript="Rscript"
fi
timestamp() { date '+%Y-%m-%d %H:%M:%S'; }
echo "$(timestamp) | job_id=${SLURM_JOB_ID:-local} task_id=${SLURM_ARRAY_TASK_ID:-local} state=START"
runtime_dir=$(mktemp -d "${SLURM_TMPDIR:-/tmp}/rw_pmm_${SLURM_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID:-local}_XXXXXX")
trap 'rm -rf "$runtime_dir"' EXIT
cd "$runtime_dir"
"$rscript" "$script_dir/simulate_batch.R" \
  "$TASK_FILE" "$SLURM_ARRAY_TASK_ID" "$RAW_DIR"
echo "$(timestamp) | job_id=${SLURM_JOB_ID:-local} task_id=${SLURM_ARRAY_TASK_ID:-local} state=COMPLETE"
