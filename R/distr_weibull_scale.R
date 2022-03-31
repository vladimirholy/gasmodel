
# WEIBULL DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_weibull_scale_parameters <- function(n) {
  group_of_par_names <- c("scale", "shape")
  par_names <- c("scale", "shape")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_weibull_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_density <- suppressWarnings(stats::dweibull(y, scale = s, shape = b))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_weibull_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dweibull(y, scale = s, shape = b, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_weibull_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_mean <- s * gamma(1 + 1 / b)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_weibull_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_var <- s^2 * (gamma(1 + 2 / b) - (gamma(1 + 1 / b))^2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_weibull_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (b * (y / s)^b - b) / s
  res_score[, 2] <- -(b * (y / s)^b * log(y / s) + b * log(s / y) - 1) / b
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_weibull_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- b^2 / s^2
  res_fisher[, 1, 2] <- -(digamma(1) + 1) / s
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (digamma(1)^2 + trigamma(1) + 2 * digamma(1) + 1) / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_weibull_scale_random <- function(t, f) {
  s <- f[1]
  b <- f[2]
  res_random <- suppressWarnings(stats::rweibull(t, scale = s, shape = b))
  res_random <- matrix(res_random, nrow = t, ncol = 1)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_weibull_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  b <- 1
  s <- y_mean
  for (i in 1:1e3) {
    b <- b + (y_var - s^2 * (gamma(1 + 2 / b) - (gamma(1 + 1 / b))^2)) / ((2 * gamma(1 / b)^2 * digamma(1 + 1 / b) - 4 * b * gamma(2 / b) * digamma((2 + b) / b)) / b^4 * s^2)
    s <- y_mean / gamma(1 + 1 / b)
  }
  b <- max(b, 1e-6)
  s <- max(s, 1e-6)
  res_start <- c(s, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


