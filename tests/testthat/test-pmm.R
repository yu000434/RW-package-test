test_that("top-k probabilities are invariant to a common intercept shift", {
  x_obs <- cbind(`(Intercept)` = 1, x = seq(-2, 2, length.out = 12))
  x_mis <- cbind(`(Intercept)` = 1, x = c(-1.5, 0, 1.5))
  beta_hat <- c(0.5, 1)
  beta_dot <- c(0.7, 0.9)
  covariance <- diag(c(0.1, 0.2))
  for (k in c(1L, 5L, nrow(x_obs))) {
    matching <- topk_probability(x_obs, x_mis, beta_hat, beta_dot, covariance, k)
    shifted <- topk_probability(x_obs, x_mis, beta_hat + c(2, 0),
                                beta_dot + c(2, 0), covariance, k)
    expect_equal(matching$probability, shifted$probability)
    expect_equal(matching$derivative, shifted$derivative)
    expect_true(all(matching$derivative[["(Intercept)"]] == 0))
    expect_equal(colSums(matching$probability), rep(k, nrow(x_mis)))
    expect_true(all(vapply(matching$derivative, function(x) {
      max(abs(colSums(x))) < 1e-8
    }, logical(1))))
  }
})

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

  pooled <- pool_rw(with_rw(recorded, {
    expect_setequal(ls(all.names = TRUE), names(recorded$data))
    lm(y ~ x)
  }))
  expect_true(all(is.finite(vcov(pooled))))
  expect_true(all(diag(vcov(pooled)) > 0))
  expect_error(pool_rw(with_rw(recorded, lm(I(2 * y) ~ x))), "untransformed PMM variable")
  expect_error(pool_rw(with_rw(recorded, lm(log(abs(y) + 1) ~ x))), "untransformed PMM variable")
})

test_that("pmmrw defaults and option restrictions are explicit", {
  set.seed(25)
  x <- matrix(rnorm(60), ncol = 1, dimnames = list(NULL, "x"))
  y <- x[, 1] + rnorm(60)
  ry <- seq_along(y) <= 40
  args <- list(y = y, ry = ry, x = x, task = "train")
  for (donors in list(5L, NULL, 0L, 100L)) {
    model <- new.env()
    set.seed(26)
    implicit <- do.call(mice.impute.pmmrw, c(args, list(model = model, donors = donors)))
    explicit_model <- new.env()
    set.seed(26)
    explicit <- do.call(mice.impute.pmmrw,
                        c(args, list(wy = !ry, model = explicit_model, donors = donors)))
    k <- if (is.null(donors)) round(sum(ry) / 600 + 7) else donors
    k <- max(1L, min(k, sum(ry)))
    expect_identical(implicit, explicit)
    expect_identical(model$donor_id, explicit_model$donor_id)
    expect_equal(model$setup$donors, k)
    set.seed(26)
    expect_identical(implicit, mice::mice.impute.pmm(y, ry, x, donors = k))
  }
  for (option in list(list(exclude = 0), list(use.matcher = TRUE),
                      list(mlocal = 2L), list(matchtype = 0L), list(matchtype = 2L))) {
    expect_error(do.call(mice.impute.pmmrw, c(args, list(model = new.env()), option)),
                 "pmmrw.*requires")
  }
  expect_error(mice.impute.pmmrw(y, ry, x), 'tasks = "train"', fixed = TRUE)
  expect_error(mice.impute.pmmrw(factor(y), ry, x, task = "train", model = new.env()),
               "numeric variables only")
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

  sigma <- imp$models$a[[1L]]$sigma.dot
  tail_score <- pmm_score_binomial(fit, "a", 1L, is.na(data$a),
                                   -c(10, 11, 12) * sigma, 0, 12L)
  for (score in tail_score) {
    expect_equal(as.vector(score), numeric(length(score)), tolerance = 1e-12)
  }
})
