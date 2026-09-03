make_rw_data <- function(seed, n = 2000, type = "robins_1") {
  set.seed(seed)
  A <- rbinom(n, 1, 1 / 3)
  X <- ifelse(A == 1, runif(n, 0.8, 2.0), runif(n, 0.1, 0.8))

  eta <- switch(type,
    robins_2_1 = c(1, 1),
    robins_2_2 = c(-1, -1),
    robins_2_3 = c(2, 1),
    c(0, 0)
  )
  mu <- X + if (type == "robins_3_2") 0.5 * X^2 else 0
  sd <- X^(ifelse(A == 1, eta[2], eta[1]) / 2)
  Z <- rnorm(n, mu, sd)

  R <- rep(1L, n)
  toddler <- which(A == 1)
  if (type == "robins_3_1") {
    R[toddler] <- rbinom(length(toddler), 1, plogis(2 - 0.5 * Z[toddler]))
  } else {
    R[toddler] <- rbinom(length(toddler), 1, 0.4)
  }

  dat <- data.frame(ID = seq_len(n), Z = Z, X = X, A = A, ImputedZ = as.integer(R == 0))
  dat$Z[R == 0] <- NA
  list(dat = dat, beta0 = 1, scenario = type)
}
