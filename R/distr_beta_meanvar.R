
# BETA DISTRIBUTION / MEAN-VARIANCE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_beta_meanvar_parameters <- function(n) {
  group_of_par_names <- c("mean", "var")
  par_names <- c("mean", "var")
  par_support <- c("probability", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_beta_meanvar_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s >= m * (1 - m)] <- NA_real_
  s[s <= m * (1 - m)] <- NA_real_
  res_density <- be_silent(stats::dbeta(y, shape1 = m * (m * (1 - m) / s - 1), shape2 = (1 - m) * (m * (1 - m) / s - 1)))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_beta_meanvar_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s >= m * (1 - m)] <- NA_real_
  s[s <= m * (1 - m)] <- NA_real_
  res_loglik <- be_silent(stats::dbeta(y, shape1 = m * (m * (1 - m) / s - 1), shape2 = (1 - m) * (m * (1 - m) / s - 1), log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_beta_meanvar_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s >= m * (1 - m)] <- NA_real_
  s[s <= m * (1 - m)] <- NA_real_
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_beta_meanvar_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s >= m * (1 - m)] <- NA_real_
  s[s <= m * (1 - m)] <- NA_real_
  res_var <- s
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_beta_meanvar_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s >= m * (1 - m)] <- NA_real_
  s[s <= m * (1 - m)] <- NA_real_
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (m^2 - m + s) / ((m - 1) * s) * (digamma(m * (1 - m) / s - 1) - digamma(m * (m * (1 - m) / s - 1)) + log(y))
  res_score[, 2] <- (s^2 * (3 * m^2 - 2 * m + s)) / (m * (m - 1) * (m^2 - m + s)) * (digamma(m * (1 - m) / s - 1) - digamma((1 - m) * (m * (1 - m) / s - 1)) + log(1 - y))
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_beta_meanvar_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s >= m * (1 - m)] <- NA_real_
  s[s <= m * (1 - m)] <- NA_real_
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- (m^2 - m + s)^2 / ((m - 1)^2 * s^2) * (trigamma(m * ((m - m^2) / s - 1)) - trigamma((m - m^2) / s - 1))
  res_fisher[, 1, 2] <- (s * ((2 - 3 * m) * m - s)) / (m * ((m - 2) * m + 1)) * trigamma((m - m^2) / s - 1)
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (s^2 * (3 * m^2 - 2 * m + s)^2) / ((m - 1)^4 * m^2) * (trigamma((1 - m) * ((m - m^2) / s - 1)) - trigamma((m - m^2) / s - 1))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_beta_meanvar_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  if (s < m * (1 - m)) {
    res_random <- be_silent(stats::dbeta(t, shape1 = m * (m * (1 - m) / s - 1), shape2 = (1 - m) * (m * (1 - m) / s - 1)))
  } else {
    res_random <- rep(NA_real_, times = t)
  }
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_beta_meanvar_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  m <- max(min(y_mean, 1 - 1e-6), 1e-6)
  s <- max(min(y_var, y_mean * (1 - y_mean) - 1e-6), 1e-6)
  res_start <- c(m, s)
  return(res_start)
}
# ------------------------------------------------------------------------------


