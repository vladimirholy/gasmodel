
# LOMAX DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_lomax_scale_parameters <- function(n) {
  group_of_par_names <- c("scale", "shape")
  par_names <- c("scale", "shape")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_lomax_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_density <- b / s / (1 + y / s)^(b + 1)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_lomax_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_loglik <- log(b) - log(s) - (b + 1) * log(1 + y / s)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_lomax_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_mean <- s / (b - 1)
  res_mean[b <= 1] <- NA_real_
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_lomax_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_var <- s^2 * b / (b - 1)^2 / (b - 2)
  res_var[b <= 2] <- NA_real_
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_lomax_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t , ncol = 2L)
  res_score[, 1] <- (b * y / s - 1) / s / (y / s + 1)
  res_score[, 2] <- 1 / b - log(y / s + 1)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_lomax_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- b / s^2 / (b + 2)
  res_fisher[, 1, 2] <- -1 / s / (b + 1)
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- 1 / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_lomax_scale_random <- function(t, f) {
  s <- f[1]
  b <- f[2]
  res_random <- be_silent(s * ((1 - stats::runif(t))^(-1 / b) - 1))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_lomax_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  s <- max(y_mean * (y_var + y_mean^2) / (y_var - y_mean^2), 1e-6)
  b <- max(2 * y_var / (y_var - y_mean^2), 1e-6)
  res_start <- c(s, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


