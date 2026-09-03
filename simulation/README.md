# RW simulation

This directory contains the parametric and PMM simulations. Parametric
imputation uses the MICE `norm` and `logreg` methods with the ordinary RW
variance. PMM uses standard MICE completed values and realized donor IDs;
smoothed matching probabilities are used only for its variance score.

## Structure

```
R/                       self-contained RW and PMM method code
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
| `analysis.R` | Constructs the analysis components `U` and `tau` for `lm` and `glm`. |
| `parametric.R` | Constructs `S_mis` and `d` for the MICE `norm` and `logreg` models. |
| `pmm.R` | Constructs the PMM imputation components, matching derivatives, PMM `kappa`, and donor information. |
| `pmmrw.R` | Runs ordinary MICE PMM and records the realized donor IDs and fitted PMM quantities. |
| `fit.R` | Evaluates the analysis and imputation components in every completed dataset. |
| `variance.R` | Assembles the ordinary RW variance or its PMM donor-source version. |
| `results.R` | Extracts variance components, Rubin-rule results, and donor-reuse summaries. |

The component flow is:

```text
analysis.R                 -> U, tau
parametric.R or pmm.R      -> S_mis, d
pmm.R                      -> PMM kappa and donor information when needed
fit.R                      -> components for all completed datasets
variance.R                 -> Omega, kappa, alpha, d_bar, and the RW variance
```

The two imputation files share the same `S_mis` and `d` interface. A simulation
may use both: in the GS PMM setting, `A` uses `pmm.R`, while the logistic
imputation of `D` uses `parametric.R`.

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

The method code depends on `mice` for imputation but does not require the old
`rw` package.

## Requirements

The simulations use `mvtnorm` and the official `amices/mice` development
version that records fitted imputation models:

```
repository: https://github.com/amices/mice
commit: c9b67ae1cd54784267a01f1b70d93a40d509a5de
version: 3.17.3.9000
```

The fixed commit is required because the simulations use `tasks = "train"` to
retain the fitted `norm`, `logreg`, and PMM quantities. No code or result is
read from `MI_4_14` or the old `rw` package.
