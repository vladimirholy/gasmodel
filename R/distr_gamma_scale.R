
# GAMMA DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_gamma_scale_parameters <- function(n) {
  group_of_par_names <- c("scale", "shape")
  par_names <- c("scale", "shape")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_gamma_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_density <- suppressWarnings(stats::dgamma(y, scale = s, shape = a))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_gamma_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dgamma(y, scale = s, shape = a, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_gamma_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_mean <- a * s
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_gamma_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_var <- a * s^2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_gamma_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (y - a * s) / s^2
  res_score[, 2] <- log(y / s) - digamma(a)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_gamma_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- a / s^2
  res_fisher[, 1, 2] <- 1 / s
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- trigamma(a)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_gamma_scale_random <- function(t, f) {
  s <- f[1]
  a <- f[2]
  res_random <- suppressWarnings(stats::rgamma(t, scale = s, shape = a))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_gamma_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  s <- max(y_var / y_mean, 1e-6)
  a <- max(y_mean^2 / y_var, 1e-6)
  res_start <- c(s, a)
  return(res_start)
}
# ------------------------------------------------------------------------------


