# rw

`rw` computes Robins-Wang variance estimates for analyses of multiply imputed
data. The first package version supports MICE `norm` and `logreg` imputation
with `lm` and binomial `glm` analysis models.

```r
imp <- mice::mice(data, method = method, tasks = "train")
fit <- with_rw(imp, lm(y ~ x))
pool_rw(fit)
```

The PMM variance method remains under development and is not included in the
package interface.
