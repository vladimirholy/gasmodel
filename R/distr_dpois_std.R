
# DOUBLE POISSON DISTRIBUTION / STANDARD PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_dpois_std_parameters <- function(n) {
  group_of_par_names <- c("mean", "dispersion")
  par_names <- c("mean", "dispersion")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_dpois_std_density <- function(y, f) {
  res_loglik <- distr_dpois_std_loglik(y = y, f = f)
  res_density <- exp(res_loglik)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_dpois_std_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  c <- 1 + (1 - s) / (12 * s * m) * (1 + 1 / (s * m))
  res_loglik <- log(s) / 2 - s * m - log(c)
  res_loglik[y > 0] <- res_loglik[y > 0] - y + s * y + y * log(y) - lfactorial(y) + y * s * log(m / y)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_dpois_std_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_dpois_std_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_var <- m / s
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_dpois_std_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- s / m * (y - m)
  res_score[, 2] <- 1 / (2 * s) - m
  res_score[y > 0, 2] <- res_score[y > 0, 2] + y * (1 + log(m) - log(y))
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_dpois_std_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- s / m
  res_fisher[, 2, 2] <- 1 / (2 * s^2)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_dpois_std_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  c <- 1 + (1 - s) / (12 * s * m) * (1 + 1 / (s * m))
  y <- 0:max(1000, 10 * m, 100 * sqrt(m / s))
  l <- log(s) / 2 - s * m - log(c) - y + s * y
  l[-1] <- l[-1] + y[-1] * log(y[-1]) - lfactorial(y[-1]) + y[-1] * s * log(m / y[-1])
  p <- exp(l)
  res_random <- sample(y, size = t, replace = TRUE, prob = p)
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_dpois_std_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  m <- max(y_mean, 1e-6)
  s <- max(y_mean / y_var, 1e-6)
  res_start <- c(m, s)
  return(res_start)
}
# ------------------------------------------------------------------------------


