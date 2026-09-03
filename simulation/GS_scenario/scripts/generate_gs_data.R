# Generates one GS simulation dataset and its true target coefficient.

make_gs_data <- function(seed, obs = 4000, subsample_n = 1000, threshold = 2) {
  set.seed(seed)

  X <- mvtnorm::rmvnorm(obs, c(0, 0), matrix(c(1, -0.25, -0.25, 1), 2, 2))
  X1 <- X[, 1]
  X2 <- X[, 2]
  Xmat <- cbind(1, X1, X2)

  A.star <- rnorm(obs, Xmat %*% c(1, 0, 0), 1)
  D.star <- rbinom(obs, 1, plogis(Xmat %*% c(-3, 0, 0) + 0.5 * A.star))
  A <- rnorm(obs, Xmat %*% c(0, -1, 0.5) + 0.9 * A.star + 0.5 * D.star, 1)
  D <- rbinom(obs, 1, plogis(Xmat %*% c(-5.5, -2, 1) + 5 * D.star + 0.5 * A))

  dat_full <- data.frame(
    ID = seq_len(obs),
    X1 = X1,
    X2 = X2,
    A.star = A.star,
    D.star = D.star,
    A = A,
    D = D,
    Intercept = 1
  )

  validated <- sample(obs, subsample_n)
  dat <- dat_full
  dat$A[-validated] <- NA
  dat$D[-validated] <- NA

  beta0 <- coef(glm(D ~ A, data = subset(dat_full, A > threshold), family = binomial()))[["A"]]
  list(dat = dat, beta0 = beta0)
}
