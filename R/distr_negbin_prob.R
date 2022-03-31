
# NEGATIVE BINOMIAL DISTRIBUTION / PROBABILISTIC PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_negbin_prob_parameters <- function(n) {
  group_of_par_names <- c("probability", "size")
  par_names <- c("probability", "size")
  par_support <- c("probability", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_negbin_prob_density <- function(y, f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  n <- f[, 2, drop = FALSE]
  res_density <- suppressWarnings(stats::dnbinom(y, prob = p, size = n))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function -------------------------------------------------------------
distr_negbin_prob_loglik <- function(y, f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  n <- f[, 2, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dnbinom(y, prob = p, size = n, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_negbin_prob_mean <- function(f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  n <- f[, 2, drop = FALSE]
  res_mean <- n * (1 - p) / p
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_negbin_prob_var <- function(f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  n <- f[, 2, drop = FALSE]
  res_var <- n * (1 - p) / p^2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_negbin_prob_score <- function(y, f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  n <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (n * p + p * y - n) / (p^2 - p)
  res_score[, 2] <- log(p) + digamma(n + y) - digamma(n)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_negbin_prob_fisher <- function(f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  n <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- -n / (p^3 - p^2)
  res_fisher[, 1, 2] <- -1 / p
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- trigamma(n) - trigamma(n + n * (1 - p) / p)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_negbin_prob_random <- function(t, f) {
  p <- f[1]
  n <- f[2]
  res_random <- suppressWarnings(stats::rnbinom(t, prob = p, size = n))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_negbin_prob_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  p <- max(min(y_mean / y_var, 1 - 1e-6), 1e-6)
  n <- max(y_mean^2 / (y_var - y_mean), 1e-6)
  res_start <- c(p, n)
  return(res_start)
}
# ------------------------------------------------------------------------------


