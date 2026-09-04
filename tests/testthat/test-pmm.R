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
