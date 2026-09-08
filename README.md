
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
remotes::install_github("amices/mice@dev")
remotes::install_github("yu000434/RW-package-test")
```

## Parametric imputation

The first step is to impute your data using `mice`. For parametric
imputation, the package supports `method = "norm"` and
`method = "logreg"`. When calling `mice`, set `tasks = "train"` to save
the fitted imputation models needed by `with_rw()` and `pool_rw()`. The
example below uses the `nhanes` data from the `mice` package. After
standardizing the variables, it creates five imputed datasets using the
`norm` method.

``` r
library(rw)
library(mice)

nhanes_scaled <- as.data.frame(scale(nhanes))

set.seed(1)
imp_norm <- mice(nhanes_scaled, method = "norm", m = 5,
                 tasks = "train", print = FALSE)
```

Use `with_rw()` to fit a linear regression of `bmi` on `age` and `hyp`
in each completed dataset.

``` r
fit_norm <- with_rw(imp_norm, lm(bmi ~ age + hyp))
```

Use `pool_rw()` to combine the coefficient estimates and compute
Robins-Wang standard errors and confidence intervals.

``` r
pool_rw(fit_norm)
#> 
#> Robins-Wang pooled results
#> -------------------------
#> Number of imputations: 5
#> Sample size: 25
#> 
#>         term estimate std.error statistic p.value conf.low conf.high
#>  (Intercept) -0.04751    0.3343   -0.1421  0.8883  -0.7407    0.6457
#>          age -0.45865    0.5336   -0.8596  0.3993  -1.5652    0.6479
#>          hyp  0.17550    1.3955    0.1258  0.9011  -2.7187    3.0697
```

For comparison, Rubin’s rules can be applied to the same imputations.

``` r
pool(with(imp_norm, lm(bmi ~ age + hyp))) |>
  summary()
#>          term    estimate std.error  statistic        df   p.value
#> 1 (Intercept) -0.04751093 0.2077719 -0.2286686 15.029658 0.8222087
#> 2         age -0.45864607 0.3618038 -1.2676651  4.179919 0.2709358
#> 3         hyp  0.17549617 0.3683495  0.4764393  3.619564 0.6610792
```

## Predictive mean matching

For predictive mean matching, use `method = "pmmrw"` in `mice` and set
`tasks = "train"`. This method makes the same donor draws as ordinary
MICE PMM and saves the donor IDs and fitted matching models needed for
variance estimation. The correction accounts for repeated use of the
same donor; it does not change the imputed values. The example below
imputes `bmi` using `age` as a predictor, with the default pool of five
donors.

``` r
nhanes_pmm <- nhanes_scaled[c("age", "bmi")]
method <- make.method(nhanes_pmm)
method[] <- ""
method["bmi"] <- "pmmrw"

set.seed(2)
imp_pmm <- mice(nhanes_pmm, method = method, m = 5,
                tasks = "train", print = FALSE)
```

Next, use `with_rw()` to regress `bmi` on `age` in each completed
dataset. Because the PMM variable is the response of a linear model,
`pool_rw()` computes the PMM cross term and applies the donor-source
correction automatically.

``` r
fit_pmm <- with_rw(imp_pmm, lm(bmi ~ age))
pool_rw(fit_pmm)
#> 
#> Robins-Wang pooled results
#> -------------------------
#> Number of imputations: 5
#> Sample size: 25
#> 
#>         term estimate std.error statistic p.value conf.low conf.high
#>  (Intercept)   0.1042    0.2422    0.4301  0.6711  -0.3969   0.60524
#>          age  -0.3036    0.1944   -1.5618  0.1320  -0.7057   0.09851
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

The package includes `giganti_data`, a simulated dataset based on
[Giganti and Shepherd (2020)](https://doi.org/10.1093/aje/kwaa153). The
example below uses all 4,000 observations. Validated values of `A` and
`D` are available for 1,000 observations and are missing for the
remainder. The analysis fits a logistic regression of `D` on `A` among
observations with `A > 2`.

``` r
data("giganti_data", package = "rw")
gs <- giganti_data
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
#> 
#> Robins-Wang pooled results
#> -------------------------
#> Number of imputations: 2
#> Sample size: 4000
#> 
#>         term estimate std.error statistic   p.value conf.low conf.high
#>  (Intercept)  -2.9450   0.25162   -11.704 1.211e-31  -3.4382   -2.4519
#>            A   0.6941   0.08243     8.421 3.728e-17   0.5326    0.8557
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
#> 
#> Robins-Wang pooled results
#> -------------------------
#> Number of imputations: 2
#> Sample size: 4000
#> 
#>         term estimate std.error statistic   p.value conf.low conf.high
#>  (Intercept)  -3.3101    0.3081   -10.742 6.449e-27  -3.9141   -2.7062
#>            A   0.8014    0.1001     8.008 1.166e-15   0.6053    0.9975
```

## Current scope

Parametric imputation supports `norm` and `logreg` with unweighted
`lm()` or binomial `glm()` analyses using the logit link.

PMM supports one numeric imputed variable with fully observed matching
predictors. Use `pool_rw()` for an untransformed PMM response in `lm()`,
or `pmm_kappa_binomial()` for the binomial setting shown above. Other
PMM analyses require an expected-score function passed to `pmm_kappa()`.
See `?pool_rw` and `?pmm_kappa` for the supported model structures, and
`?mice.impute.pmmrw` for PMM options.
