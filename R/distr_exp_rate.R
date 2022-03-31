
# EXPONENTIAL DISTRIBUTION / RATE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_exp_rate_parameters <- function(n) {
  group_of_par_names <- c("rate")
  par_names <- c("rate")
  par_support <- c("positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_exp_rate_density <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  res_density <- suppressWarnings(stats::dexp(y, rate = r))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_exp_rate_loglik <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dexp(y, rate = r, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_exp_rate_mean <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  res_mean <- 1 / r
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_exp_rate_var <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  res_var <- 1 / r^2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_exp_rate_score <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 1L)
  res_score[, 1] <- 1 / r - y
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_exp_rate_fisher <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 1L, 1L))
  res_fisher[, 1, 1] <- 1 / r^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_exp_rate_random <- function(t, f) {
  r <- f[1]
  res_random <- suppressWarnings(stats::rexp(t, rate = r))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_exp_rate_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  r <- max(1 / y_mean, 1e-6)
  res_start <- r
  return(res_start)
}
# ------------------------------------------------------------------------------


