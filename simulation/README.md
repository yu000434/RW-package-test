# RW simulation

This directory contains the parametric and PMM simulations. Parametric
imputation uses the MICE `norm` and `logreg` methods with the ordinary RW
variance. PMM uses standard MICE completed values and realized donor IDs;
smoothed matching probabilities are used only for its variance score.

Both simulations load the package source directly from `../R/`. There is no
separate copy of the method code and no installed `rw` package is required.

## Structure

```
../R/                    shared package core
results.R                result extraction and diagnostics
GS_scenario/scripts/     GS generator and simulation scripts
RW_scenario/scripts/     RW generator and simulation scripts
GS_scenario/tasks_*.csv  fixed GS task tables
RW_scenario/tasks_*.csv  fixed RW task tables
results/<run-id>/        raw results and summaries
```

GS and RW have separate data generators and one-replicate runners because
their data-generating and analysis models differ. Both use the same files in
`../R/` and write the same result columns. Copy the whole `RW` project to run
on another machine, not just this `simulation/` directory.

The existing `rb_*` raw-result names are retained so the frozen PMM files
remain compatible. In parametric rows, these columns contain the ordinary RW
variance and standard error.

## Calculation

`with_rw()` fits each completed dataset. PMM cross terms use `pmm_kappa()`
for RW or `pmm_kappa_binomial()` for GS. `pool_rw()` computes the final
variance in both the parametric and PMM paths.

`results.R` saves Rubin-rule results, variance components, score norms,
donor reuse and matching errors. It reconstructs variance components only
for reporting and checks that their sum equals the package variance.
Matching diagnostics call the package's probability code; they do not define
a second PMM method. This adds one probability calculation per imputation to
retain the diagnostic fields without changing the package interface. The
original task tables, seeds and result columns are unchanged.

## Target coefficient

`beta0` is the population coefficient of the complete-data analysis model.
For the GS generator and `A > 2`, it is approximately `0.885377976781648`,
obtained by numerical integration of the population logistic score.
The same value is used for all sample sizes and imputation settings.
For RW, `beta0 = 1`, except for S3b: the quadratic mean gives the
no-intercept linear slope `1.689707792207792`. S3a and S3b are extensions
of the original RW settings.

Earlier saved results used a per-sample GS coefficient and `1` for RW S3b.
Their `beta0`, bias and true-target coverage are not comparable to new runs.
Empirical-centered coverage remains a separate variance diagnostic.

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
