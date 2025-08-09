
# GENERALIZED ERROR DISTRIBUTION / MEAN-SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_ged_meanscale_parameters <- function(n) {
  group_of_par_names <- c("mean", "scale", "shape")
  par_names <- c("mean", "scale", "shape")
  par_support <- c("real", "positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_ged_meanscale_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_density <- be_silent(b / 2 / s / gamma(1 / b) * exp(-(abs(y - m) / s)^b))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_ged_meanscale_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_loglik <- be_silent(log(b) - log(2) - log(s) - lgamma(1 / b) - (abs(y - m) / s)^b)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_ged_meanscale_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_ged_meanscale_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_var <- s^2 * gamma(3 / b) / gamma(1 / b)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_ged_meanscale_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- (y != m) * b * (1 / s)^b * (y - m) * abs(y - m)^(b - 2)
  res_score[, 2] <- b * (1 / s)^(b + 1) * abs(y - m)^b - 1 / s
  res_score[, 3] <- 1 / b + s^(-b) * abs(y - m)^b * (log(s) - log(abs(y - m)) + 1e-6) + digamma(1 / b) / b^2
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_ged_meanscale_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  for (i in 1:t) {
    y_range <- m[i, ] + c(-3, 3) * s[i, ]^2 * gamma(3 / b[i, ]) / gamma(1 / b[i, ])
    y_seq <- seq(from = y_range[1], to = y_range[2], length.out = 1e2)
    y_by <- y_seq[2] - y_seq[1]
    y_seq <- as.matrix(y_seq[abs(y_seq - m[1, ]) > 1e-9])
    for (j in 1:length(y_seq)) {
      density <- as.vector(distr_ged_meanscale_density(y_seq[j, , drop = FALSE], f[i, , drop = FALSE]))
      score <- distr_ged_meanscale_score(y_seq[j, , drop = FALSE], f[i, , drop = FALSE])
      res_fisher[i, , ] <- res_fisher[i, , ] + t(score) %*% score * density * y_by
    }
  }
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_ged_meanscale_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  b <- f[3]
  res_random <- be_silent(m + sample(c(-1, 1), replace = TRUE, size = t) * (stats::rgamma(t, scale = s^b, shape = 1 / b))^(1 / b))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_ged_meanscale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_kurt <- mean((y - y_mean)^4) / (mean((y - y_mean)^2)^2) - 3
  m <- y_mean
  b <- 5 * (y_kurt < -1) + 2 * (y_kurt >= -1 & y_kurt < 1) + 1 * (y_kurt >= 1 & y_kurt <= 5) + 0.7 * (y_kurt > 5)
  s <- sqrt(y_var * gamma(1 / b) / gamma(3 / b))
  res_start <- c(m, s, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


