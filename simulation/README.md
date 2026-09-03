# RW PMM simulation

This is the simulation workspace for `MICE_soft_PMM_donor_source_research`.
It is separate from the future `rw` package. MICE PMM completed values and
recorded donor IDs are retained exactly. Soft top-k probabilities are used
only to construct the donor-source-corrected variance estimate.

## Structure

```
R/                       common PMM variance code
GS_scenario/PMM/         GS task table and scripts
RW_scenario/PMM/         RW task tables and scripts
results/<run-id>/        raw files, logs, summary, and provenance
reports/                  paper-style results report
tests/                    MICE parity and deterministic smoke checks
```

GS and RW use separate data generators and one-replicate runners because
their data-generating and analysis models differ. Both use the same files in
`R/` and write the same result columns.

## Common R files

| File | Purpose |
|---|---|
| `pmmrw.R` | Calls the ordinary MICE PMM matcher, returns the same copied donor values, and records realized donor IDs and fitted PMM quantities. |
| `topk.R` | Constructs soft top-k donor probabilities and their derivatives with respect to the PMM coefficients. |
| `rb.R` | Computes the Rao--Blackwellized PMM cross term. It contains the RW conditional score and the GS integrated logistic score. |
| `rw_variance.R` | Applies donor-source clustering and assembles the corrected RW sandwich variance. |
| `metrics.R` | Extracts variance components, Rubin-rule quantities, and donor-reuse summaries for the simulation output. |

These files contain the research estimator. They do not generate GS or RW
data and do not submit jobs.

## Scenario scripts

Each `PMM/scripts/` directory follows the same workflow:

| File | Purpose |
|---|---|
| `generate_*_data.R` | Generates one simulated dataset and its true target value. |
| `run_one.R` | Runs one seed: generate data, impute with MICE PMM, fit the analysis, and calculate RB and Rubin-rule results. |
| `make_tasks.R` | Creates a task CSV for a requested simulation design. |
| `audit_tasks.R` | Checks task IDs, parameter values, replicate counts, and non-overlapping seeds before submission. |
| `simulate_batch.R` | Runs all replicates assigned to one task, prints replicate progress, and atomically writes one raw CSV. |
| `run_array_task.sh` | Slurm array entry point that passes one array task to `simulate_batch.R`. |
| `submit_cell.sh` or `submit_run.sh` | Finds missing task IDs, freezes source hashes, and submits only those tasks. |
| `audit_results.R` | Checks that all expected tasks, seeds, scenarios, and result rows exist after the jobs finish. |
| `combine_results.R` | Combines fully audited raw files and calculates coverage summaries and Monte Carlo intervals. |

GS uses `submit_cell.sh` because one submission selects one `(m, k)` cell.
RW uses `submit_run.sh` because each replicate evaluates all six RW scenarios
for one fixed `(n, m, k)` combination.

The active task tables are stored directly in the corresponding `PMM/`
directory. `GS_scenario/PMM/tasks_n2000.csv` is one table containing all 16
GS `(m, k)` cells. RW currently has two tables:
`tasks_n1000_m50_k5.csv` and `tasks_n150_m20_k5.csv`. A task table is
immutable once a run starts. Its run ID must state the relevant parameters,
such as `gs_n2000_m25_k5_v1` or `rw_n1000_m50_k5_v1`.

Normally, a user calls a submit script, then `audit_results.R`, and finally
`combine_results.R`. `run_array_task.sh` and `simulate_batch.R` are internal
execution steps called by Slurm.

## Runs

For a 2,500-replicate run, use 50 tasks with 50 replicates per task. Audit the
task table before submission. After all tasks are terminal, audit the raw
files before writing `summary.csv`. Do not combine partial output.

```sh
Rscript GS_scenario/PMM/scripts/audit_tasks.R GS_scenario/PMM/tasks_n2000.csv
Rscript GS_scenario/PMM/scripts/audit_results.R TASK_FILE RAW_DIR CELL_ID
Rscript GS_scenario/PMM/scripts/combine_results.R TASK_FILE RAW_DIR SUMMARY_FILE CELL_ID
```

Use `DRY_RUN=1` with a submit script before a new Slurm submission. Submit
scripts create `logs/` and `provenance/` next to the supplied `RAW_DIR`, so
use a centralized path such as `results/gs_n2000_m25_k5_v1/raw`.

## Checks

Run these after any source cleanup:

```sh
Rscript tests/test_mice_pmm_parity.R
Rscript tests/smoke.R
```

`smoke.R` compares one GS and one RW replicate against frozen raw results.
These checks verify implementation parity; they are not new simulation
evidence.
