
<!-- README.md is generated from README.Rmd. Please edit README.Rmd. -->

# `rw`: Robins-Wang variance estimation for multiple imputation

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

## Current scope

The package supports `norm` and `logreg` imputation with `lm()` or
binomial `glm()` analysis. Automatic PMM pooling currently supports one
numeric `pmmrw` variable used as the response of `lm()`. More complex
PMM analyses can provide their cross-term matrix through the `pmm_kappa`
argument to `pool_rw()`.
