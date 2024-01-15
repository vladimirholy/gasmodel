
# GENERALIZED GAMMA DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_burr_scale_parameters <- function(n) {
  group_of_par_names <- c("scale", "shape1", "shape2")
  par_names <- c("scale", "shape1", "shape2")
  par_support <- c("positive", "positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_burr_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_density <- a * b / s * (y / s)^(a - 1) / (1 + (y / s)^a)^(b + 1)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_burr_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_loglik <- log(a) + log(b) - log(s) + (a - 1) * (log(y) - log(s)) - (b + 1) * log(1 + (y / s)^a)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_burr_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_mean <- s * b * beta(b - 1 / a, 1 + 1 / a)
  res_mean[a <= 1] <- NA_real_
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_burr_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_var <- s^2 * b * beta(b - 2 / a, 1 + 2 / a) - s^2 * b^2 * beta(b - 1 / a, 1 + 1 / a)^2
  res_var[a <= 2] <- NA_real_
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_burr_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t , ncol = 3L)
  res_score[, 1] <- a * (b * (y / s)^a - 1) / s / ((y / s)^a + 1)
  res_score[, 2] <- 1 / a - (b * (y / s)^a - 1) * log(y / s) / ((y / s)^a + 1)
  res_score[, 3] <- 1 / b - log((y / s)^a + 1)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_burr_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- a^2 * b / s^2 / (b + 2)
  res_fisher[, 1, 2] <- -b * (1 + digamma(1) - digamma(b + 1)) / s / (b + 2)
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 1, 3] <- -a / s / (b + 1)
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 2, 2] <- (1 + b / (b + 2) * (pi^2 / 6 + digamma(1)^2 + 2 * digamma(1) - 2 * (digamma(1) + 1) * digamma(b + 1) + digamma(b + 1)^2 + trigamma(b + 1))) / a^2
  res_fisher[, 2, 3] <- (1 + digamma(1) - digamma(b)) / a / (b + 1)
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <- 1 / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_burr_scale_random <- function(t, f) {
  s <- f[1]
  a <- f[2]
  b <- f[3]
  res_random <- be_silent(s * ((1 - stats::runif(t))^(-1 / b) - 1)^(1 / a))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_burr_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_square <- mean(y^2, na.rm = TRUE)
  a <- 3
  s <- y_mean * sqrt(3) * 3 / 2 / pi
  for (i in 1:1e3) {
    a <- a + (y_square - s^2 * 2 * pi / a / sin(2 * pi / a)) / (2 * pi * s^2 * (2 * pi * pracma::cot(2 * pi / a) - a) * pracma::csc(2 * pi / a) / a^3)
    s <- y_mean * sin(pi / a) / pi * a
  }
  a <- max(a, 1e-6)
  s <- max(s, 1e-6)
  b <- 1
  res_start <- c(s, a, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


