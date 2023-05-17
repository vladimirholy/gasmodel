
# BETA DISTRIBUTION / CONCENTRATION PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_beta_conc_parameters <- function(n) {
  group_of_par_names <- c("conc", "conc")
  par_names <- c("conc1", "conc2")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_beta_conc_density <- function(y, f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_density <- be_silent(stats::dbeta(y, shape1 = a, shape2 = b))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_beta_conc_loglik <- function(y, f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_loglik <- be_silent(stats::dbeta(y, shape1 = a, shape2 = b, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_beta_conc_mean <- function(f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_mean <- a / (a + b)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_beta_conc_var <- function(f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_var <- (a * b) / ((a + b)^2 * (a + b + 1))
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_beta_conc_score <- function(y, f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- digamma(a + b) - digamma(a) + log(y)
  res_score[, 2] <- digamma(a + b) - digamma(b) + log(1 - y)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_beta_conc_fisher <- function(f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- trigamma(a) - trigamma(a + b)
  res_fisher[, 1, 2] <- -trigamma(a + b)
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- trigamma(b) - trigamma(a + b)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_beta_conc_random <- function(t, f) {
  a <- f[1]
  b <- f[2]
  res_random <- be_silent(stats::dbeta(t, shape1 = a, shape2 = b))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_beta_conc_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  a <- max(y_mean * (y_mean * (1 - y_mean) / y_var - 1), 1e-6)
  b <- max((1 - y_mean) * (y_mean * (1 - y_mean) / y_var - 1), 1e-6)
  res_start <- c(a, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


