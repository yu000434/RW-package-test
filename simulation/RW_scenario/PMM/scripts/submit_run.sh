#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: submit_run.sh TASK_FILE RAW_DIR" >&2
  exit 2
fi

task_file=$(cd "$(dirname "$1")" && pwd)/$(basename "$1")
raw_dir=$2
script_dir=$(cd "$(dirname "$0")" && pwd)
pmm_dir=$(cd "$script_dir/.." && pwd)
run_dir=$(dirname "$raw_dir")
log_dir="$run_dir/logs"
provenance_dir="$run_dir/provenance"
walltime="${WALLTIME:-01:00:00}"
qos_args=()
if [[ -n "${QOS:-}" ]]; then
  qos_args=(--qos="$QOS")
fi

task_ids=$(awk -F, 'NR > 1 {print $1}' "$task_file" | paste -sd, -)
if [[ -z "$task_ids" ]]; then
  echo "No RW tasks found." >&2
  exit 1
fi

missing_ids=()
IFS=, read -r -a ids <<< "$task_ids"
for task_id in "${ids[@]}"; do
  printf -v result_name 'rw_task%04d.csv' "$task_id"
  [[ -f "$raw_dir/$result_name" ]] || missing_ids+=("$task_id")
done
if [[ ${#missing_ids[@]} -eq 0 ]]; then
  echo "All RW task outputs already exist."
  exit 0
fi

array_ids=$(IFS=,; echo "${missing_ids[*]}")
if [[ "${DRY_RUN:-0}" == 1 ]]; then
  echo "DRY_RUN array_ids=$array_ids raw_dir=$raw_dir"
  exit 0
fi

mkdir -p "$raw_dir" "$log_dir" "$provenance_dir"
raw_dir=$(cd "$raw_dir" && pwd)
current_hashes=$(mktemp)
sha256sum "$task_file" "$script_dir/simulate_batch.R" "$script_dir/run_array_task.sh" \
  "$script_dir/run_one.R" "$script_dir/generate_rw_data.R" \
  "$pmm_dir/../../R/"*.R > "$current_hashes"
if [[ ! -e "$provenance_dir/task_table.csv" ]]; then
  cp "$task_file" "$provenance_dir/task_table.csv"
  mv "$current_hashes" "$provenance_dir/source_sha256.txt"
elif ! cmp -s "$task_file" "$provenance_dir/task_table.csv"; then
  rm -f "$current_hashes"
  echo "The task table differs from the frozen provenance copy." >&2
  exit 1
elif ! cmp -s "$current_hashes" "$provenance_dir/source_sha256.txt"; then
  rm -f "$current_hashes"
  echo "The source files differ from the frozen provenance hash." >&2
  exit 1
else
  rm -f "$current_hashes"
fi

submission=$(sbatch --job-name="rw_pmm" --array="$array_ids" \
  --cpus-per-task=1 --mem=2G --time="$walltime" \
  "${qos_args[@]}" \
  --output="$log_dir/%A_%a.log" \
  --export="ALL,TASK_FILE=$task_file,RAW_DIR=$raw_dir" \
  "$script_dir/run_array_task.sh")
printf '%s | array_ids=%s walltime=%s | %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" \
  "$array_ids" "$walltime" "$submission" | tee -a "$provenance_dir/submission.log"
