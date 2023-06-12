
# GENERALIZED GAMMA DISTRIBUTION / RATE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_gengamma_rate_parameters <- function(n) {
  group_of_par_names <- c("rate", "shape1", "shape2")
  par_names <- c("rate", "shape1", "shape2")
  par_support <- c("positive", "positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_gengamma_rate_density <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_density <- b * r * (y * r)^(a * b - 1) * exp(-(y * r)^b) / gamma(a)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_gengamma_rate_loglik <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_loglik <- log(b) + log(r) + (a * b - 1) * log(y * r) - (y * r)^b - lgamma(a)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_gengamma_rate_mean <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_mean <- 1 / r * gamma(a + 1 / b) / gamma(a)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_gengamma_rate_var <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_var <- 1 / r^2 * (gamma(a + 2 / b) / gamma(a) - (gamma(a + 1 / b) / gamma(a))^2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_gengamma_rate_score <- function(y, f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- -((r * y)^b * b - a * b) / r
  res_score[, 2] <- b * log(r * y) - digamma(a)
  res_score[, 3] <- -((r * y)^b * b * log(r * y) - a * b * log(r * y) - 1) / b
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_gengamma_rate_fisher <- function(f) {
  t <- nrow(f)
  r <- f[, 1, drop = FALSE]
  a <- f[, 2, drop = FALSE]
  b <- f[, 3, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- a * b^2 / r^2
  res_fisher[, 1, 2] <- -b / r
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 1, 3] <- (a * digamma(a) + 1) / r
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 2, 2] <- trigamma(a)
  res_fisher[, 2, 3] <- -digamma(a) / b
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <- (a * digamma(a)^2 + a * trigamma(a) + 2 * digamma(a) + 1) / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_gengamma_rate_random <- function(t, f) {
  r <- f[1]
  a <- f[2]
  b <- f[3]
  res_random <- 1 / r * be_silent(stats::rgamma(t, shape = a, rate = 1))^(1 / b)
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_gengamma_rate_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_cubic <- mean(y^3, na.rm = TRUE)
  y_quartic <- mean(y^4, na.rm = TRUE)
  b <- 1
  a <- y_mean^2 / y_var
  r <- y_mean / y_var
  for (i in 1:1e3) {
    b <- b + (y_cubic - 1 / r^3 * gamma(a + 3 / b) / gamma(a)) / (-3 * gamma(a + 3 / b) * digamma(a + 3 / b) / b^2 / r^3 / gamma(a))
    a <- a + (y_quartic - 1 / r^4 * (gamma(a + 4 / b) / gamma(a))) / (gamma(a + 4 / b) * (digamma(a + 4 / b) - digamma(a)) / r^4 / gamma(a))
    r <- 1 / y_mean * gamma(a + 1 / b) / gamma(a)
  }
  b <- max(b, 1e-6)
  a <- max(a, 1e-6)
  r <- max(r, 1e-6)
  res_start <- c(r, a, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


