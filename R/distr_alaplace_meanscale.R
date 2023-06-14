
# ASYMMETRIC LAPLACE DISTRIBUTION / MEAN-SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_alaplace_meanscale_parameters <- function(n) {
  group_of_par_names <- c("mean", "scale", "asym")
  par_names <- c("mean", "scale", "asym")
  par_support <- c("real", "positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_alaplace_meanscale_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  a <- f[, 3, drop = FALSE]
  res_density <- 1 / s / (a + 1 / a) * exp(-abs(y - m) / s * a^sign(y - m))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_alaplace_meanscale_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  a <- f[, 3, drop = FALSE]
  res_loglik <- -log(s) - log(a + 1 / a) - abs(y - m) / s * a^sign(y - m)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_alaplace_meanscale_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  a <- f[, 3, drop = FALSE]
  res_mean <- m + s * (1 / a - a)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_alaplace_meanscale_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  a <- f[, 3, drop = FALSE]
  res_var <- s^2 * (1 / a^2 + a^2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_alaplace_meanscale_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  a <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- sign(y - m) * a^sign(y - m) / s
  res_score[, 2] <- abs(y - m) * a^sign(y - m) / s^2 - 1 / s
  res_score[, 3] <- -(y - m) * a^(sign(y - m) - 1) / s + (1 - a^2) / (a * (1 + a^2))
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_alaplace_meanscale_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  a <- f[, 3, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- 1 / s^2
  res_fisher[, 1, 3] <- -2 / s / (1 + a^2)
  res_fisher[, 2, 2] <- 1 / s^2
  res_fisher[, 2, 3] <- -1 / s / a * (1 - a^2) / (1 + a^2)
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <- 1 / a^2 + 4 / (1 + a^2) ^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_alaplace_meanscale_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  a <- f[3]
  u <- stats::runif(n = t, min = -a, max = 1 / a)
  res_random <- m - s * sign(u) / a^sign(u) * log(1 - abs(u) * a^sign(u))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_alaplace_meanscale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_skew <- mean((y - y_mean)^3) / y_var^(3 / 2)
  a <- 1
  try_and_be_silent(a <- stats::uniroot(function (x) { 2 * (1 - x^6) / (1 + x^4)^(3 / 2) - min(max(y_skew, -1.99), 1.99)}, lower = 0.01, upper = 100, tol = 1e-3)$root)
  s <- max(sqrt(y_var / (1 + a^4) * a^2), 1e-6)
  m <- y_mean - s * (1 - a^2) / a
  res_start <- c(m, s, a)
  return(res_start)
}
# ------------------------------------------------------------------------------


