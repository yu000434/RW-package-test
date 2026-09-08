# Generate the Giganti and Shepherd (2020) example data.
# Adapted from LucyMcGowan/rw: data-raw/giganti-data.R.
# Run from the package root with mvtnorm installed.

set.seed(455)
n <- 4000L
X <- mvtnorm::rmvnorm(n, mean = c(0, 0),
                      sigma = matrix(c(1, -0.25, -0.25, 1), 2, 2))
Xmat <- cbind(1, X)
beta_star <- c(1, 0, 0)
gamma_star <- c(-3, 0, 0, 0.5)
beta <- c(0, -1, 0.5, 0.9, 0.5)
gamma <- c(-5.5, -2, 1, 0, 5, 0.5)

# Error-prone measurements, available for every observation.
A.star <- rnorm(n, Xmat %*% beta_star, 1)
eta <- Xmat %*% gamma_star[1:3] + A.star * gamma_star[4]
D.star <- rbinom(n, 1, exp(eta) / (1 + exp(eta)))

# Validated measurements, retained only in the sampled 1,000 observations.
eta <- Xmat %*% beta[1:3] + A.star * beta[4] + D.star * beta[5]
A <- rnorm(n, eta, 1)
eta <- Xmat %*% gamma[1:3] + A.star * gamma[4] +
  D.star * gamma[5] + A * gamma[6]
D <- rbinom(n, 1, exp(eta) / (1 + exp(eta)))
sampled <- seq_len(n) %in% sample(seq_len(n), 1000L, replace = FALSE)

giganti_data <- data.frame(ID = seq_len(n), X1 = X[, 1], X2 = X[, 2],
                           A.star, D.star, A, D, Intercept = 1)
giganti_data$A <- ifelse(sampled, A, NA)
giganti_data$D <- factor(ifelse(sampled, D, NA), levels = c(0, 1))
save(giganti_data, file = "data/giganti_data.rda", compress = "xz", version = 2)
