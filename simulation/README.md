# RW PMM simulation

This directory reproduces the simulations for
`MICE_soft_PMM_donor_source_research`. MICE PMM completed values and realized
donor IDs are retained exactly. Soft top-k probabilities are used only in the
donor-source-corrected variance estimate.

## Structure

```
R/                       common PMM variance code
GS_scenario/PMM/         GS scripts and fixed task table
RW_scenario/PMM/         RW scripts and fixed task tables
results/<run-id>/        raw results and summaries
reports/                  paper-style report
tests/                    parity and smoke tests
```

GS and RW have separate data generators and one-replicate runners because
their data-generating and analysis models differ. Both use the same files in
`R/` and write the same result columns.

## Method files

| File | Purpose |
|---|---|
| `pmmrw.R` | Runs ordinary MICE PMM and records the realized donor IDs and fitted PMM quantities. |
| `soft_pmm.R` | Calculates soft-PMM donor probabilities and their derivatives. |
| `pmm_score.R` | Calculates the PMM cross term for the RW and GS analyses. |
| `variance.R` | Applies donor-source clustering and assembles the corrected RW variance. |
| `results.R` | Extracts variance components, Rubin-rule results, and donor-reuse summaries. |

## Simulation scripts

Each scenario has four scripts:

| File | Purpose |
|---|---|
| `generate_*_data.R` | Generates one simulated dataset and its true target value. |
| `run_one.R` | Runs one seed from data generation through imputation, analysis, and variance estimation. |
| `run.R` | Runs one row of a fixed task table and writes one raw CSV. |
| `summarize.R` | Checks that the requested raw results are complete and calculates the paper summary and Monte Carlo intervals. |

The fixed GS table, `GS_scenario/PMM/tasks_n2000.csv`, contains all 16
`(m, k)` cells. The two fixed RW tables correspond to `n = 1000, m = 50,
k = 5` and `n = 150, m = 20, k = 5`.

Run different task IDs in parallel with any scheduler. The simulation code
does not depend on Slurm:

```sh
Rscript GS_scenario/PMM/scripts/run.R TASK_FILE TASK_ID RAW_DIR
Rscript GS_scenario/PMM/scripts/summarize.R TASK_FILE RAW_DIR SUMMARY_FILE CELL_ID

Rscript RW_scenario/PMM/scripts/run.R TASK_FILE TASK_ID RAW_DIR
Rscript RW_scenario/PMM/scripts/summarize.R TASK_FILE RAW_DIR SUMMARY_FILE
```

## Checks

```sh
Rscript tests/test_mice_pmm_parity.R
Rscript tests/smoke.R
```

The parity test verifies that `pmmrw` preserves ordinary MICE completed
values and records the correct donor IDs. The smoke test compares one GS and
one RW replicate against frozen raw results.
