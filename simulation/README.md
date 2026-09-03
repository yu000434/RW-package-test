# RW simulation

This directory contains the parametric and PMM simulations. Parametric
imputation uses the MICE `norm` and `logreg` methods with the ordinary RW
variance. PMM uses standard MICE completed values and realized donor IDs;
smoothed matching probabilities are used only for its variance score.

## Structure

```
R/                       common variance and PMM code
GS_scenario/scripts/     GS generator and simulation scripts
RW_scenario/scripts/     RW generator and simulation scripts
GS_scenario/tasks_*.csv  fixed GS task tables
RW_scenario/tasks_*.csv  fixed RW task tables
results/<run-id>/        raw results and summaries
reports/                  paper-style report
tests/                    parity and smoke tests
```

GS and RW have separate data generators and one-replicate runners because
their data-generating and analysis models differ. Both use the same files in
`R/` and write the same result columns.

The existing `rb_*` raw-result names are retained so the frozen PMM files
remain compatible. In parametric rows, these columns contain the ordinary RW
variance and standard error.

## Method files

| File | Purpose |
|---|---|
| `pmmrw.R` | Runs ordinary MICE PMM and records the realized donor IDs and fitted PMM quantities. |
| `pmm_score.R` | Calculates PMM matching probabilities, their derivatives, and the PMM cross term for the RW and GS analyses. |
| `variance.R` | Assembles the ordinary RW variance or its PMM donor-source version. |
| `results.R` | Extracts variance components, Rubin-rule results, and donor-reuse summaries. |

## Simulation scripts

Each scenario has four scripts:

| File | Purpose |
|---|---|
| `generate_*_data.R` | Generates one simulated dataset and its true target value. |
| `run_one.R` | Runs one seed from data generation through imputation, analysis, and variance estimation. |
| `run.R` | Runs one row of a fixed task table and writes one raw CSV. |
| `summarize.R` | Checks that the requested raw results are complete and calculates the paper summary and Monte Carlo intervals. |

`tasks_pmm.csv` contains the 16 GS `(m, k)` cells at `n = 2000`.
`tasks_parametric.csv` contains the certified GS parametric grid at
`n = 2000, 4000, 6000, 8000` and `m = 5, 25, 50, 100`. The RW task tables
contain the two settings `n = 150, m = 20` and `n = 1000, m = 50`; PMM uses
`k = 5`. Every full cell has 2,500 replications.

Run different task IDs in parallel with any scheduler. The simulation code
does not depend on Slurm:

```sh
Rscript GS_scenario/scripts/run.R TASK_FILE TASK_ID RAW_DIR
Rscript GS_scenario/scripts/summarize.R TASK_FILE RAW_DIR SUMMARY_FILE CELL_ID

Rscript RW_scenario/scripts/run.R TASK_FILE TASK_ID RAW_DIR
Rscript RW_scenario/scripts/summarize.R TASK_FILE RAW_DIR SUMMARY_FILE
```

## Checks

```sh
Rscript tests/test_mice_pmm_parity.R
Rscript tests/smoke.R
```

The parity test verifies that `pmmrw` preserves ordinary MICE completed
values and records the correct donor IDs. The smoke test checks PMM against
the frozen current results and parametric MICE against the certified 2026-07-18
results.
