test_that("pmmrw preserves MICE values and records donors", {
  set.seed(21)
  data <- data.frame(x = rnorm(120))
  data$y <- 1 + 2 * data$x + rnorm(120)
  data$y[sample(120, 40)] <- NA
  method <- mice::make.method(data)
  method[] <- ""
  method["y"] <- "pmm"
  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0L
  pred["y", "x"] <- 1L

  set.seed(22)
  ordinary <- mice::mice(data, m = 3, method = method, predictorMatrix = pred,
                         blots = list(y = list(donors = 5)), print = FALSE)
  method["y"] <- "pmmrw"
  set.seed(22)
  recorded <- mice::mice(data, m = 3, method = method, predictorMatrix = pred,
                         blots = list(y = list(donors = 5)), tasks = "train",
                         print = FALSE)

  donor_id <- extract_donor_id(recorded, "y")
  missing <- which(is.na(data$y))
  for (p in seq_len(recorded$m)) {
    standard <- mice::complete(ordinary, p)$y
    completed <- mice::complete(recorded, p)$y
    expect_identical(completed, standard)
    expect_equal(completed[missing], data$y[donor_id[missing, p]])
  }

  pooled <- pool_rw(with_rw(recorded, lm(y ~ x)))
  expect_true(all(is.finite(vcov(pooled))))
  expect_true(all(diag(vcov(pooled)) > 0))
})

test_that("PMM works with downstream logreg and binomial analysis", {
  set.seed(23)
  data <- data.frame(x1 = rnorm(240), x2 = rnorm(240))
  data$a <- 0.5 + data$x1 - 0.5 * data$x2 + rnorm(240)
  data$d <- factor(rbinom(240, 1, plogis(-1 + data$a + 0.5 * data$x1)))
  missing <- sample(240, 120)
  data$a[missing] <- NA
  data$d[missing] <- NA
  method <- mice::make.method(data)
  method[] <- ""
  method[c("a", "d")] <- c("pmmrw", "logreg")
  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0L
  pred["a", c("x1", "x2")] <- 1L
  pred["d", c("x1", "x2", "a")] <- 1L

  set.seed(24)
  imp <- mice::mice(data, m = 2, method = method, predictorMatrix = pred,
                    tasks = "train", print = FALSE)
  fit <- with_rw(imp, glm(d ~ a, subset = a > 0, family = binomial()))
  kappa <- pmm_kappa_binomial(fit, "a", threshold = 0, quadrature_order = 12)
  pooled <- pool_rw(fit, pmm_kappa = kappa)

  expect_equal(dim(kappa), c(2, 3))
  expect_true(all(is.finite(vcov(pooled))))
  expect_true(all(diag(vcov(pooled)) > 0))
})
