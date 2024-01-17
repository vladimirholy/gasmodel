
# BIRNBAUM-SAUNDERS DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_bisa_scale_parameters <- function(n) {
  group_of_par_names <- c("scale", "shape")
  par_names <- c("scale", "shape")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_bisa_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_density <- sqrt(s / y) * (1 + s / y) * exp((2 - y / s - s / y) / (2 * a^2)) / (2 * a * s * sqrt(2 * pi))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_bisa_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_loglik <- log(s / y) / 2 + log(1 + s / y) + ((2 - y / s - s / y) / (2 * a^2)) - log(2 * a * s * sqrt(2 * pi))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_bisa_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_mean <- s * (1 + a^2 / 2)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_bisa_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_var <- (s * a)^2 * (1 + 5 * a^2 / 4)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_bisa_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- y / (2 * a^2 * s^2) - 1 / (2 * a^2 * y) + 1 / (s + y) - 1 / (2 * s)
  res_score[, 2] <- y / (a^3 * s) + s / (a^3 * y) - (2 + a^2) / a^3
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_bisa_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- (1 + a * (a * sqrt(pi / 2) - pi * exp(2 / a^2) * (1 - stats::pnorm(2 / a))) / sqrt(2 * pi)) / (a * s)^2
  res_fisher[, 2, 2] <- 2 / a^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_bisa_scale_random <- function(t, f) {
  s <- f[1]
  a <- f[2]
  x <- be_silent(stats::rnorm(n = t, mean = 0, sd = a / 2))
  res_random <- s * (1 + 2 * x^2 + 2 * x * sqrt((1 + x^2)))
  res_random <- matrix(res_random, nrow = t, ncol = 1)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_bisa_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_med <- stats::median(y, na.rm = TRUE)
  s <- max(y_med, 1e-6)
  a <- sqrt(max(2 * (y_mean / y_med - 1), 1e-12))
  res_start <- c(s, a)
  return(res_start)
}
# ------------------------------------------------------------------------------


