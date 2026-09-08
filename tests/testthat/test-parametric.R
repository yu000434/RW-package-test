test_that("normal imputation works with linear analysis", {
  set.seed(11)
  data <- data.frame(x = rnorm(80))
  data$y <- 1 + 2 * data$x + rnorm(80)
  data$y[sample(80, 24)] <- NA
  method <- mice::make.method(data)
  method[] <- ""
  method["y"] <- "norm"
  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0L
  pred["y", "x"] <- 1L
  imp <- mice::mice(data, m = 2, method = method, predictorMatrix = pred,
                    tasks = "train", print = FALSE)

  fit <- with_rw(imp, {
    expect_setequal(ls(all.names = TRUE), names(imp$data))
    lm(y ~ x)
  })
  pooled <- pool_rw(fit)

  expect_s3_class(fit, "rw_fit")
  expect_s3_class(pooled, "rw_pool")
  expect_named(coef(pooled), c("(Intercept)", "x"))
  expect_equal(dim(vcov(pooled)), c(2, 2))
  expect_true(all(is.finite(vcov(pooled))))

  output <- capture.output(printed <- withVisible(print(pooled)))
  expect_true("Number of imputations: 2" %in% output)
  expect_true("Sample size: 80" %in% output)
  expect_identical(printed$value, pooled)
  expect_false(printed$visible)
})

test_that("normal and logistic imputation work with binomial analysis", {
  set.seed(12)
  data <- data.frame(x = rnorm(200))
  data$a <- 0.5 * data$x + rnorm(200)
  data$d <- factor(rbinom(200, 1, plogis(-0.5 + 0.8 * data$a + 0.2 * data$x)))
  missing <- sample(200, 60)
  data$a[missing] <- NA
  data$d[missing] <- NA
  method <- mice::make.method(data)
  method[] <- ""
  method[c("a", "d")] <- c("norm", "logreg")
  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0L
  pred["a", "x"] <- 1L
  pred["d", c("x", "a")] <- 1L
  imp <- mice::mice(data, m = 2, method = method, predictorMatrix = pred,
                    tasks = "train", print = FALSE)

  pooled <- pool_rw(with_rw(imp, glm(d ~ a, family = binomial())))

  expect_named(coef(pooled), c("(Intercept)", "a"))
  expect_true(all(is.finite(vcov(pooled))))
  expect_s3_class(summary(pooled), "data.frame")
})
