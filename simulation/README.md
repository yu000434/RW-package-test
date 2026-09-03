# RW PMM simulation

This is the simulation workspace for `MICE_soft_PMM_donor_source_research`.
It is separate from the future `rw` package. MICE PMM completed values and
recorded donor IDs are retained exactly; the soft top-k calculation is used
only in the donor-source-corrected variance estimate.

```
R/method/             active PMM variance code
GS_scenario/PMM/      GS data, runner, scripts, and task tables
RW_scenario/PMM/      RW data, runner, scripts, and task tables
results/<run-id>/     raw task files, logs, summary.csv, and provenance
reports/              paper-style results report
tests/                MICE parity and deterministic smoke checks
```

GS and RW keep separate data generators and one-replicate runners because
their data-generating and analysis models differ. They use the same PMM
variance code and output columns.

## Runs

A task table is immutable once a run starts. Its run ID must include the
relevant parameters, for example `gs_n2000_m25_k5_v1` or
`rw_n1000_m50_k5_v1`. One GS submission selects one `(m, k)` cell. One RW
submission evaluates all six RW scenarios for a fixed `(n, m, k)`.

For a 2,500-replicate run, use 50 tasks with 50 replicates per task. Audit a
task table before submission. After all tasks are terminal, audit the raw
files and then write `summary.csv`. Do not combine partial output.

```
Rscript GS_scenario/PMM/scripts/audit_tasks.R TASK_FILE
Rscript GS_scenario/PMM/scripts/audit_results.R TASK_FILE RAW_DIR CELL_ID
Rscript GS_scenario/PMM/scripts/combine_results.R TASK_FILE RAW_DIR SUMMARY_FILE CELL_ID
```

Use `DRY_RUN=1` with a submit script before a new Slurm submission. The
submit scripts create logs and provenance next to the supplied `RAW_DIR`, so
use a centralized path such as `results/gs_n2000_m25_k5_v1/raw`.

## Checks

Run these after any source cleanup:

```
Rscript tests/test_mice_pmm_parity.R
Rscript tests/smoke.R
```

`smoke.R` compares one GS and one RW replicate against frozen raw results.
These checks are implementation checks, not new simulation evidence.
