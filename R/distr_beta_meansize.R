
# BETA DISTRIBUTION / MEAN-SIZE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_beta_meansize_parameters <- function(n) {
  group_of_par_names <- c("mean", "size")
  par_names <- c("mean", "size")
  par_support <- c("probability", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_beta_meansize_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_density <- be_silent(stats::dbeta(y, shape1 = m * v, shape2 = (1 - m) * v))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_beta_meansize_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_loglik <- be_silent(stats::dbeta(y, shape1 = m * v, shape2 = (1 - m) * v, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_beta_meansize_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_beta_meansize_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_var <- m * (1 - m) / (v + 1)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_beta_meansize_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- v / (1 - m) * (digamma(v) - digamma(m * v) + log(y))
  res_score[, 2] <- digamma(v) - digamma(v - v * m) + log(1 - y)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_beta_meansize_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- v^2 / (1 - m)^2 * (trigamma(m * v) - trigamma(v))
  res_fisher[, 1, 2] <- v / (1 - m) * (-trigamma(v))
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- trigamma(v - m * v) - trigamma(v)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_beta_meansize_random <- function(t, f) {
  m <- f[1]
  v <- f[2]
  res_random <- be_silent(stats::dbeta(t, shape1 = m * v, shape2 = (1 - m) * v))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_beta_meansize_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  m <- max(min(y_mean, 1 - 1e-6), 1e-6)
  v <- max(y_mean * (1 - y_mean) / y_var - 1, 1e-6)
  res_start <- c(m, v)
  return(res_start)
}
# ------------------------------------------------------------------------------


