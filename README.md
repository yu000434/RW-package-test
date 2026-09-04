# rw: Robins-Wang variance estimation for multiple imputation

`rw` computes Robins-Wang variance estimates for analyses of multiply imputed
data. The estimator accounts for uncertainty from both the imputation and
analysis models and provides an alternative to Rubin's rules.

The methodology is based on [Robins and Wang
(2000)](https://doi.org/10.1093/biomet/87.1.113).

## Installation

The package requires the development version of `mice` that records fitted
imputation models through `tasks = "train"`.

```r
remotes::install_github("amices/mice@c9b67ae1cd54784267a01f1b70d93a40d509a5de")
remotes::install_github("yu000434/rw-package-test")
```

## Example

This example imputes three incomplete variables using normal working models
and fits a linear regression in each completed dataset.

```r
library(mice)
library(rw)

data <- as.data.frame(scale(mice::nhanes))

method <- make.method(data)
method[] <- ""
method[c("bmi", "hyp", "chl")] <- "norm"

set.seed(1)
imp <- mice(
  data,
  m = 5,
  method = method,
  tasks = "train",
  print = FALSE
)

fit <- with_rw(imp, lm(bmi ~ age + hyp))
pooled <- pool_rw(fit)
summary(pooled)
```

The pooled coefficients and variance-covariance matrix are also available
directly:

```r
coef(pooled)
vcov(pooled)
```

## Current scope

The current package interface supports:

- MICE `norm` and `logreg` imputation models;
- linear analysis models fitted with `lm()`;
- binomial analysis models fitted with `glm()`;
- model subsets specified in the analysis expression.

Predictive mean matching variance estimation is still being validated and is
not included in this package version.
