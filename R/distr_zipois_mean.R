
# ZERO-INFLATED POISSON DISTRIBUTION / MEAN PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_zipois_mean_parameters <- function(n) {
  group_of_par_names <- c("rate", "inflation")
  par_names <- c("rate", "inflation")
  par_support <- c("positive", "probability")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_zipois_mean_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_density <- (y == 0L) * p + (1 - p) * suppressWarnings(stats::dpois(y, lambda = m))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function -------------------------------------------------------------
distr_zipois_mean_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_loglik <- matrix(0, nrow = t, ncol = 1)
  res_loglik[y > 0L, ] <- log(1 - p[y > 0L, ]) + suppressWarnings(stats::dpois(y[y > 0L, ], lambda = m[y > 0L], log = TRUE))
  res_loglik[y == 0L, ] <- log(p[y == 0L, ] + (1 - p[y == 0L, ]) * suppressWarnings(stats::dpois(0L, lambda = m[y == 0L, ])))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_zipois_mean_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_mean <- (1 - p) * m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_zipois_mean_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_var <- m * (1 - p) * (1 + p * m)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_zipois_mean_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (y == 0L) * (p - 1) / (p * (exp(m) - 1) + 1) + (y > 0L) * (y - m) / m
  res_score[, 2] <- (y == 0L) * (exp(m) - 1)/(p * (exp(m) - 1) + 1) + (y > 0L) * 1 / (p - 1)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_zipois_mean_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  p <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- ((m - exp(m) + 1) * p^2 - (m - exp(m) + 2) * p + 1) / ((m * exp(m) - m) * p + m)
  res_fisher[, 1, 2] <- -1 / (p * (exp(m) - 1) + 1)
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- -(exp(m) - 1) / (p^2 * (exp(m) - 1) - p * (exp(m) - 2) - 1)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_zipois_mean_random <- function(t, f) {
  m <- f[1]
  p <- f[2]
  res_random <- sample(c(0L, NA_real_), size = t, replace = TRUE, prob = c(p, 1 - p))
  res_random[is.na(res_random)] <- suppressWarnings(stats::rpois(sum(is.na(res_random)), lambda = m))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_zipois_mean_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_zero <- mean(y == 0L, na.rm = TRUE)
  p <- 0
  m <- y_mean
  for (i in 1:1e3) {
    p <- (y_zero - exp(-m)) / (1 - exp(-m))
    m <- y_mean / (1 - p)
  }
  p <- max(min(p, 1 - 1e-6), 1e-6)
  m <- max(m, 1e-6)
  res_start <- c(m, p)
  return(res_start)
}
# ------------------------------------------------------------------------------


