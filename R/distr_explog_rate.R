
# EXPONENTIAL-LOGARITHMIC DISTRIBUTION / RATE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_explog_rate_parameters <- function(n) {
  group_of_par_names <- c("rate", "shape")
  par_names <- c("rate", "shape")
  par_support <- c("positive", "probability")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_explog_rate_density <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_density <- r / log(p) / (1 - exp(r * y) / (1 - p))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_explog_rate_loglik <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_loglik <- log(r) - log(-log(p)) - log(exp(r * y) / (1 - p) - 1)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_explog_rate_mean <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_mean <- -copula::polylog(1 - p, 2) / r / log(p)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_explog_rate_var <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_var <- -2 * copula::polylog(1 - p, 3, n.sum = 1000) / r^2 / log(p) - (copula::polylog(1 - p, 2) / r / log(p))^2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_explog_rate_score <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- 1 / r - y - ((1 - p) * y * exp(-r * y)) / (1 - (1 - p) * exp(-r * y))
  res_score[, 2] <- -1 / (p * log(p)) - 1 / (1 - p) - exp(-r * y) / (1 - (1 - p) * exp(-r * y))
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_explog_rate_fisher <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- -copula::polylog(1 - p, 2) / r^2 / log(p)
  res_fisher[, 1, 2] <- (1 - p + p * log(p)) / (2 * r * p * (1 - p) * log(p))
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- -(log(p) + 1) / (p * log(p))^2 + 1 / (1 - p)^2 + ((1 - 4 * p) / p^2 - 2 * log(p) + 3) / (2 * (1 - p)^2 * log(p))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_explog_rate_random <- function(t, f) {
  r <- f[1]
  p <- f[2]
  res_random <- be_silent(log((1 - p) / (1 - p^(stats::runif(t)))) / r)
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_explog_rate_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  n <- nrow(y)
  p <- 0.5
  r <- max(-copula::polylog(1 - p, 2) / y_mean / log(p), 1e-6)
  for (i in 1:1e2) {
    p <- max(min(-n * (1 - p) / (log(p) * sum(1 / (1 - (1 - p) * exp(-r * y)))), 1 - 1e-6), 1e-6)
    r <- max(n / (sum(y / (1 - (1 - p) * exp(-r * y)))), 1e-6)
  }
  res_start <- c(r, p)
  return(res_start)
}
# ------------------------------------------------------------------------------


