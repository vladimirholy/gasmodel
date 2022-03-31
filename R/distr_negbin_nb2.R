
# NEGATIVE BINOMIAL DISTRIBUTION / NB2 PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_negbin_nb2_parameters <- function(n) {
  group_of_par_names <- c("mean", "dispersion")
  par_names <- c("mean", "dispersion")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_negbin_nb2_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_density <- suppressWarnings(stats::dnbinom(y, mu = m, size = 1 / s))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_negbin_nb2_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dnbinom(y, mu = m, size = 1 / s, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_negbin_nb2_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_negbin_nb2_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_var <- m * (1 + s * m)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_negbin_nb2_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (y - m) / m / (1 + s * m)
  res_score[, 2] <- (y - m) / s / (1 + s * m) + (log(1 + s * m) + digamma(1 / s) - digamma(y + 1 / s)) / s^2
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_negbin_nb2_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- 1 / m / (1 + s * m)
  res_fisher[, 2, 2] <- 2 * log(1 + s * m) / s^3 - m / (s^2 + s^3 * m) + 2 * digamma(1 / s) / s^3 + trigamma(1 / s) / s^4 - 2 * digamma(m + 1 / s) / s^3 - trigamma(m + 1 / s) / s^4
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_negbin_nb2_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  res_random <- suppressWarnings(stats::rnbinom(t, mu = m, size = 1 / s))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_negbin_nb2_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  m <- max(y_mean, 1e-6)
  s <- max((y_var - y_mean) / y_mean^2, 1e-6)
  res_start <- c(m, s)
  return(res_start)
}
# ------------------------------------------------------------------------------


