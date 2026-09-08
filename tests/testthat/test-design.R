test_that("factor designs match the matrices supplied to MICE", {
  norm <- getFromNamespace("mice.impute.norm", "mice")
  logreg <- getFromNamespace("mice.impute.logreg", "mice")
  local_mocked_bindings(
    mice.impute.norm = function(y, ry, x, model, ...) {
      model$test_x <- cbind(`(Intercept)` = 1, x)
      norm(y, ry, x, model = model, ...)
    },
    mice.impute.logreg = function(y, ry, x, model, ...) {
      model$test_x <- cbind(`(Intercept)` = 1, x)
      logreg(y, ry, x, model = model, ...)
    }, .package = "mice")

  set.seed(31)
  data <- data.frame(x = rnorm(240),
                     group = factor(rep(c("c", "a", "b"), 80), levels = c("c", "b", "a")),
                     binary = factor(rep(c("no", "yes"), 120)))
  contrasts(data$group) <- contr.sum(3)
  data$d <- factor(rbinom(240, 1, plogis(data$x)), labels = c("no", "yes"))
  data$y <- data$x + as.integer(data$group) + as.integer(data$d) + rnorm(240)
  data$d[1:60] <- NA
  data$y[61:120] <- NA
  method <- c(x = "", group = "", binary = "", d = "logreg", y = "norm")
  pred <- mice::make.predictorMatrix(data)
  pred[,] <- 0
  pred["d", c("x", "group", "binary")] <- 1
  pred["y", c("x", "group", "binary", "d")] <- 1
  imp <- mice::mice(data, m = 2, maxit = 2, method = method, predictorMatrix = pred,
                    visitSequence = c("d", "y"), tasks = "train", print = FALSE)
  fit <- with_rw(imp, lm(y ~ x + group + binary + d))
  for (p in seq_len(imp$m)) {
    completed <- mice::complete(imp, p)
    for (variable in c("d", "y")) {
      model <- imp$models[[variable]][[p]]
      expect_equal(unname(imputation_design(completed, model)), unname(model$test_x))
    }
    expect_identical(fit$results[[p]]$model$model$d, completed$d)
    expect_equal(coef(fit$results[[p]]$model), coef(lm(y ~ x + group + binary + d, completed)))
  }
  expect_true(all(is.finite(vcov(pool_rw(fit)))))
  logistic <- with_rw(imp, glm(d ~ x + group, family = binomial()))
  for (p in seq_len(imp$m)) {
    model <- logistic$results[[p]]$model
    y <- as.integer(model$model$d) - 1L
    expect_equal(as.vector(logistic$results[[p]]$U),
                 as.vector(model.matrix(model) * (y - fitted(model))))
  }
  expect_true(all(is.finite(vcov(pool_rw(logistic)))))
})

test_that("PMM factor predictors agree with explicit dummy columns", {
  pmmrw <- mice.impute.pmmrw
  local_mocked_bindings(mice.impute.pmmrw = function(y, ry, x, model, ...) {
    model$test_x <- cbind(`(Intercept)` = 1, x)
    pmmrw(y, ry, x, model = model, ...)
  }, .package = "rw")
  set.seed(32)
  data <- data.frame(x = rnorm(180), group = factor(rep(c("a", "b", "c"), 60)))
  data$y <- data$x + as.integer(data$group) + rnorm(180)
  data$y[1:60] <- NA
  dummy <- data.frame(x = data$x, groupb = as.integer(data$group == "b"),
                      groupc = as.integer(data$group == "c"), y = data$y)
  impute <- function(data, method) {
    methods <- mice::make.method(data)
    methods[] <- ""
    methods["y"] <- method
    pred <- mice::make.predictorMatrix(data)
    pred[,] <- 0
    pred["y", setdiff(names(data), "y")] <- 1
    mice::mice(data, method = methods, predictorMatrix = pred, m = 2, maxit = 2,
                seed = 33, tasks = if (method == "pmm") "impute" else "train", print = FALSE)
  }
  imp <- impute(data, "pmmrw")
  ordinary <- impute(data, "pmm")
  encoded <- impute(dummy, "pmmrw")
  expect_identical(imp$imp$y, ordinary$imp$y)
  expect_identical(imp$imp$y, encoded$imp$y)
  expect_identical(extract_donor_id(imp, "y"), extract_donor_id(encoded, "y"))
  for (p in seq_len(imp$m)) {
    model <- imp$models$y[[p]]
    expect_equal(unname(imputation_design(data, model)), unname(model$test_x))
  }
  fit <- with_rw(imp, lm(y ~ x + group))
  encoded_fit <- with_rw(encoded, lm(y ~ x + groupb + groupc))
  expect_equal(pmm_kappa(fit, "y", pmm_score_lm),
               pmm_kappa(encoded_fit, "y", pmm_score_lm))
  expect_equal(vcov(pool_rw(fit)), vcov(pool_rw(encoded_fit)))
})

test_that("binomial PMM handles factor predictors and labelled responses", {
  set.seed(34)
  data <- data.frame(x = rnorm(240), group = factor(rep(c("a", "b", "c"), 80)))
  data$a <- data$x + as.integer(data$group) + rnorm(240)
  data$d <- factor(rbinom(240, 1, plogis(-2 + data$a)), labels = c("no", "yes"))
  data$a[1:80] <- NA
  data$d[41:120] <- NA
  impute <- function(data) {
    methods <- mice::make.method(data)
    methods[] <- ""
    methods[c("a", "d")] <- c("pmmrw", "logreg")
    pred <- mice::make.predictorMatrix(data)
    pred[,] <- 0
    predictors <- setdiff(names(data), c("a", "d"))
    pred["a", predictors] <- 1
    pred["d", c(predictors, "a")] <- 1
    mice::mice(data, method = methods, predictorMatrix = pred, m = 2, maxit = 2,
                seed = 35, tasks = "train", print = FALSE)
  }
  encoded <- data.frame(x = data$x, groupb = as.integer(data$group == "b"),
                        groupc = as.integer(data$group == "c"), a = data$a,
                        d = factor(as.integer(data$d) - 1L, levels = c(0, 1)))
  imp <- impute(data)
  numeric_imp <- impute(encoded)
  fit <- with_rw(imp, glm(d ~ a, subset = a > 0, family = binomial()))
  numeric_fit <- with_rw(numeric_imp, glm(d ~ a, subset = a > 0, family = binomial()))
  kappa <- pmm_kappa_binomial(fit, "a", threshold = 0, quadrature_order = 12)
  numeric_kappa <- pmm_kappa_binomial(numeric_fit, "a", threshold = 0, quadrature_order = 12)
  expect_equal(kappa, numeric_kappa)
  expect_true(all(is.finite(kappa)))
  expect_equal(vcov(pool_rw(fit, kappa)), vcov(pool_rw(numeric_fit, numeric_kappa)))
})
