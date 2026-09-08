
<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

# `rw`: Robins-Wang variance estimation for multiple imputation

<!-- badges: start -->

[![R-CMD-check](https://github.com/yu000434/RW-package-test/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yu000434/RW-package-test/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `rw` package computes [Robins-Wang variance
estimates](https://doi.org/10.1093/biomet/87.1.113) for analyses of
multiply imputed data. It supports parametric chained-equation
imputation and predictive mean matching with donor-source correction.

## Installation

The package requires the development version of `mice` that records
fitted imputation models through `tasks = "train"`.

``` r
remotes::install_github("amices/mice@c9b67ae1cd54784267a01f1b70d93a40d509a5de")
remotes::install_github("yu000434/RW-package-test")
```

## Parametric imputation

This example imputes the incomplete variables in the `nhanes` data using
normal working models.

``` r
library(rw)
library(mice)

nhanes_scaled <- as.data.frame(scale(nhanes))
method <- make.method(nhanes_scaled)
method[] <- ""
method[c("bmi", "hyp", "chl")] <- "norm"

set.seed(1)
imp_norm <- mice(nhanes_scaled, method = method, m = 5,
                 tasks = "train", print = FALSE)
fit_norm <- with_rw(imp_norm, lm(bmi ~ age + hyp))
pool_rw(fit_norm)
#> Robins-Wang pooled results
#>         term    estimate std.error   statistic   p.value  conf.low conf.high
#>  (Intercept) -0.02127991 0.7101366 -0.02996595 0.9763644 -1.494013  1.451453
#>          age -0.68329820 1.5840270 -0.43136777 0.6703986 -3.968369  2.601773
#>          hyp  0.15921513 2.1031229  0.07570415 0.9403387 -4.202395  4.520825
```

For comparison, Rubin’s rules can be applied to the same imputations.

``` r
pool(with(imp_norm, lm(bmi ~ age + hyp))) |>
  summary()
#>          term    estimate std.error   statistic        df    p.value
#> 1 (Intercept) -0.02127991 0.2201983 -0.09663976 16.243158 0.92419483
#> 2         age -0.68329820 0.2886778 -2.36699283 12.046683 0.03552134
#> 3         hyp  0.15921513 0.3108677  0.51216365  5.073478 0.63004620
```

## Predictive mean matching

The `pmmrw` method uses the ordinary MICE donor draw and records the
selected donor. The completed values are unchanged. The matching
derivative is used only for variance estimation, and `pool_rw()`
accounts for repeated use of the same donor.

Use `tasks = "train"` with numeric PMM variables and `matchtype = 1`.
The options `exclude = NULL`, `use.matcher = FALSE`, and `mlocal = 1`
remain at their defaults. The donor pool size is controlled by `donors`
(default 5).

``` r
nhanes_pmm <- nhanes_scaled[c("age", "bmi")]
method <- make.method(nhanes_pmm)
method[] <- ""
method["bmi"] <- "pmmrw"

set.seed(2)
imp_pmm <- mice(nhanes_pmm, method = method, m = 5,
                tasks = "train", print = FALSE)
fit_pmm <- with_rw(imp_pmm, lm(bmi ~ age))
pool_rw(fit_pmm)
#> Robins-Wang pooled results
#>         term   estimate std.error  statistic   p.value   conf.low  conf.high
#>  (Intercept)  0.1041708 0.2422202  0.4300667 0.6711499 -0.3968998 0.60524153
#>          age -0.3035738 0.1943690 -1.5618422 0.1319808 -0.7056568 0.09850922
```

Using the same seed with ordinary MICE PMM gives the same completed
values.

``` r
method["bmi"] <- "pmm"
set.seed(2)
imp_standard <- mice(nhanes_pmm, method = method, m = 5, print = FALSE)
identical(imp_pmm$imp$bmi, imp_standard$imp$bmi)
#> [1] TRUE
```

The donor IDs used in each imputation are available when needed.

``` r
head(extract_donor_id(imp_pmm, "bmi"))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    8    8    8   19    9
#> [2,]   NA   NA   NA   NA   NA
#> [3,]    7   15    8   22   25
#> [4,]   17    2   22   25    8
#> [5,]   NA   NA   NA   NA   NA
#> [6,]   24   20    8   24    5
```

## Giganti and Shepherd example

The package includes the simulated `giganti_data`. The first 500
observations are used here to keep the example short. Both `A` and `D`
are incomplete.

``` r
data("giganti_data", package = "rw")
gs <- giganti_data[seq_len(500), ]
predictors <- c("X1", "X2", "A.star", "D.star")

pred <- make.predictorMatrix(gs)
pred[,] <- 0L
pred["A", predictors] <- 1L
pred["D", c(predictors, "A")] <- 1L
```

With parametric imputation, `A` uses a normal model and `D` uses
logistic regression.

``` r
method <- make.method(gs)
method[] <- ""
method[c("A", "D")] <- c("norm", "logreg")

set.seed(41)
imp_gs_norm <- mice(gs, m = 2, method = method, predictorMatrix = pred,
                    tasks = "train", print = FALSE)
fit_gs_norm <- with_rw(
  imp_gs_norm,
  glm(D ~ A, subset = A > 2, family = binomial())
)
pool_rw(fit_gs_norm)
#> Robins-Wang pooled results
#>         term   estimate std.error statistic      p.value   conf.low conf.high
#>  (Intercept) -3.0864612 0.6738210 -4.580536 4.637867e-06 -4.4071261 -1.765796
#>            A  0.9385706 0.2255057  4.162070 3.153760e-05  0.4965875  1.380554
```

For predictive mean matching, `pmm_kappa_binomial()` computes the PMM
cross term for the binomial analysis before the final variance is
assembled.

``` r
method[c("A", "D")] <- c("pmmrw", "logreg")

set.seed(42)
imp_gs_pmm <- mice(gs, m = 2, method = method, predictorMatrix = pred,
                   tasks = "train", print = FALSE)
fit_gs_pmm <- with_rw(
  imp_gs_pmm,
  glm(D ~ A, subset = A > 2, family = binomial())
)
kappa <- pmm_kappa_binomial(fit_gs_pmm, "A", threshold = 2)
pool_rw(fit_gs_pmm, pmm_kappa = kappa)
#> Robins-Wang pooled results
#>         term   estimate std.error statistic    p.value   conf.low  conf.high
#>  (Intercept) -2.4084915 1.1015650 -2.186427 0.02878438 -4.5675192 -0.2494637
#>            A  0.6485696 0.3549776  1.827072 0.06768899 -0.0471738  1.3443130
```

## Current scope

The package supports `norm` and `logreg` imputation with unweighted
`lm()` or binomial `glm()` analysis using the default logit link.

PMM variance estimation currently requires one numeric `pmmrw` variable
with fully observed predictors in its matching model:

- `pool_rw()` automatically handles an untransformed PMM response in
  `lm()`. Analysis predictors and any subset must depend only on fully
  observed variables.
- `pmm_kappa_binomial()` handles `glm(D ~ A, family = binomial())`,
  where `A` uses PMM and `D` uses `logreg`. The analysis can use all
  rows or the subset `A > threshold`. The `logreg` model must include
  `A` as an untransformed additive predictor; its other predictors must
  be fully observed.
- Other PMM analyses need an expected-score function supplied to
  `pmm_kappa()`. Its arguments and required output are described in
  `?pmm_kappa`.

## Source files

| File | Role |
|:---|:---|
| `R/pmmrw.R` | Standard MICE PMM draws and donor recording. |
| `R/parametric.R` | Imputation scores and influence contributions for `norm` and `logreg`. |
| `R/pmm.R` | PMM matching probabilities, derivatives, and expected-score cross terms. |
| `R/analysis.R` | Analysis-model scores and their derivatives. |
| `R/rw.R` | Fitting, RW variance assembly, pooling, and result methods. |

Function comments in `R/` are kept short; the full reference pages are
maintained in `man/`. `devtools::document()` updates `NAMESPACE` only.
Edit `README.Rmd` and run `rmarkdown::render("README.Rmd")` to update
this page and its example output. The package website is built from this
README and the reference pages.
