
# WEIBULL DISTRIBUTION / RATE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_weibull_rate_parameters <- function(n) {
  group_of_par_names <- c("rate", "shape")
  par_names <- c("rate", "shape")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_weibull_rate_density <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_density <- be_silent(stats::dweibull(y, scale = 1 / r, shape = b))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_weibull_rate_loglik <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_loglik <- be_silent(stats::dweibull(y, scale = 1 / r, shape = b, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_weibull_rate_mean <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_mean <- 1 / r * gamma(1 + 1 / b)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_weibull_rate_var <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_var <- 1 / r^2 * (gamma(1 + 2 / b) - (gamma(1 + 1 / b))^2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_weibull_rate_score <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- -((r * y)^b * b - b) / r
  res_score[, 2] <- -((r * y)^b * b * log(r * y) - b * log(r * y) - 1) / b
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_weibull_rate_fisher <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- b^2 / r^2
  res_fisher[, 1, 2] <- (digamma(1) + 1) / r
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (digamma(1)^2 + trigamma(1) + 2 * digamma(1) + 1) / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_weibull_rate_random <- function(t, f) {
  r <- f[1]
  b <- f[2]
  res_random <- be_silent(stats::rweibull(t, scale = 1 / r, shape = b))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_weibull_rate_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  b <- 1
  r <- 1 / y_mean
  for (i in 1:1e3) {
    b <- b + (y_var - 1 / r^2 * (gamma(1 + 2 / b) - (gamma(1 + 1 / b))^2)) / ((2 * gamma(1 / b)^2 * digamma(1 + 1 / b) - 4 * b * gamma(2 / b) * digamma((2 + b) / b)) / b^4 / r^2)
    r <- 1 / y_mean * gamma(1 + 1 / b)
  }
  b <- max(b, 1e-6)
  r <- max(r, 1e-6)
  res_start <- c(r, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


