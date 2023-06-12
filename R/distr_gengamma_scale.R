
# GENERALIZED GAMMA DISTRIBUTION / SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_gengamma_scale_parameters <- function(n) {
  group_of_par_names <- c("scale", "shape1", "shape2")
  par_names <- c("scale", "shape1", "shape2")
  par_support <- c("positive", "positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_gengamma_scale_density <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_density <- b / s * (y / s)^(a * b - 1) * exp(-(y / s)^b) / gamma(a)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_gengamma_scale_loglik <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_loglik <- log(b) - log(s) + (a * b - 1) * log(y / s) - (y / s)^b - lgamma(a)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_gengamma_scale_mean <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_mean <- s * gamma(a + 1 / b) / gamma(a)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_gengamma_scale_var <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_var <- s^2 * (gamma(a + 2 / b) / gamma(a) - (gamma(a + 1 / b) / gamma(a))^2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_gengamma_scale_score <- function(y, f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t , ncol = 3L)
  res_score[, 1] <- -(a * b - b * (y / s)^b) / s
  res_score[, 2] <- b * log(y / s) - digamma(a)
  res_score[, 3] <- -(a * b * log(s / y) + b * (y / s)^b * log(y / s) - 1) / b
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_gengamma_scale_fisher <- function(f) {
  t <- nrow(f)
  s <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- a * b^2 / s^2
  res_fisher[, 1, 2] <- b / s
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 1, 3] <- -(a * digamma(a) + 1) / s
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 2, 2] <- trigamma(a)
  res_fisher[, 2, 3] <- -digamma(a) / b
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <- (a * digamma(a)^2 + a * trigamma(a) + 2 * digamma(a) + 1) / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_gengamma_scale_random <- function(t, f) {
  s <- f[1]
  a <- f[2]
  b <- f[3]
  res_random <- s * be_silent(stats::rgamma(t, shape = a, scale = 1))^(1 / b)
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_gengamma_scale_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_cubic <- mean(y^3, na.rm = TRUE)
  y_quartic <- mean(y^4, na.rm = TRUE)
  b <- 1
  a <- y_mean^2 / y_var
  s <- y_var / y_mean
  for (i in 1:1e3) {
    b <- b + (y_cubic - s^3 * gamma(a + 3 / b) / gamma(a)) / (-3 * gamma(a + 3 / b) * digamma(a + 3 / b) / b^2 * s^3 / gamma(a))
    a <- a + (y_quartic - s^4 * (gamma(a + 4 / b) / gamma(a))) / (gamma(a + 4 / b) * (digamma(a + 4 / b) - digamma(a)) * s^4 / gamma(a))
    s <- y_mean / gamma(a + 1 / b) * gamma(a)
  }
  b <- max(b, 1e-6)
  a <- max(a, 1e-6)
  s <- max(s, 1e-6)
  res_start <- c(s, a, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


