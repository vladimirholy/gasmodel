
# RAYLEIGH DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_rayleigh_scale_parameters <- function(n) {
  group_of_par_names <- c("scale")
  par_names <- c("scale")
  par_support <- c("positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_rayleigh_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  res_density <- be_silent(stats::dweibull(y, scale = s * sqrt(2), shape = 2))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_rayleigh_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  res_loglik <- be_silent(stats::dweibull(y, scale = s * sqrt(2), shape = 2, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_rayleigh_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  res_mean <- s * sqrt(pi / 2)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_rayleigh_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  res_var <- s^2 * (4 - pi) / 2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_rayleigh_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 1L)
  res_score[, 1] <- (y^2 - 2 * s^2) / s^3
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_rayleigh_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 1L, 1L))
  res_fisher[, 1, 1] <- 4 / s^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_rayleigh_scale_random <- function(t, f) {
  s <- f[1]
  res_random <- be_silent(stats::rweibull(t, scale = s * sqrt(2), shape = 2))
  res_random <- matrix(res_random, nrow = t, ncol = 1)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_rayleigh_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  s <- max(y_mean * sqrt(2 / pi), 1e-6)
  res_start <- s
  return(res_start)
}
# ------------------------------------------------------------------------------


