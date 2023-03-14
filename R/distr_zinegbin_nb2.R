
# ZERO-INFLATED NEGATIVE BINOMIAL DISTRIBUTION / NB2 PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_zinegbin_nb2_parameters <- function(n) {
  group_of_par_names <- c("mean", "dispersion", "inflation")
  par_names <- c("mean", "dispersion", "inflation")
  par_support <- c("positive", "positive", "probability")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_zinegbin_nb2_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_density <- (y == 0L) * p + (1 - p) * be_silent(stats::dnbinom(y, mu = m, size = 1 / s))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_zinegbin_nb2_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  res_loglik[y > 0L, ] <- log(1 - p[y > 0L, ]) + be_silent(stats::dnbinom(y[y > 0L, ], mu = m[y > 0L, ], size = 1 / s[y > 0L, ], log = TRUE))
  res_loglik[y == 0L, ] <- log(p[y == 0L, ] + (1 - p[y == 0L, ]) * be_silent(stats::dnbinom(0L, mu = m[y == 0L, ], size = 1 / s[y == 0L, ])))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_zinegbin_nb2_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_mean <- (1 - p) * m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_zinegbin_nb2_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_var <- m * (1 - p) * (1 + p * m + s * m)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_zinegbin_nb2_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- (y == 0L) * (p - 1)*(1 / (m * s + 1))^(1 / s)/(m * p *s - ((m * p - m) * s + p - 1) * (1 / (m * s + 1))^(1 / s) + p) + (y > 0L) * (y - m) / (m^2 * s + m)
  res_score[, 2] <- (y == 0L) * (-1) * ((m * p * (log(m * s + 1) - 1) - m * (log(m * s + 1) - 1)) * s + p * log(m * s + 1) - log(m * s + 1)) * (1 / (m * s + 1))^(1 / s) / (m * p * s^3 + p * s^2 - ((m * p - m) * s^3 + (p - 1) * s^2) * (1 / (m * s + 1))^(1 / s)) + (y > 0L) * ((m * (log(m * s + 1) - 1) + m * digamma(1 / s)) * s + s * y - (m * s + 1) * digamma(y + 1 / s) + log(m * s + 1) + digamma(1 / s)) / (m * s^3 + s^2)
  res_score[, 3] <- (y == 0L) * ((1 / (m * s + 1))^(1 / s) - 1) / ((p - 1) * (1 / (m * s + 1))^(1 / s) - p) + (y > 0L) * 1 / (p - 1)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_zinegbin_nb2_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- -(p^2 + (m * p^2 - m * p) * s - ((m + 1) * p^2 - (m + 2) * p + (m * p^2 - 2 * m * p + m) * s + 1) * (1 / (m * s + 1))^(1 / s) - p)/(m^3 * p * s^2 + 2 * m^2 * p * s + m * p - ((m^3 * p - m^3) * s^2 + m * p + 2 * (m^2 * p - m^2) * s - m) * (1 / (m * s + 1))^(1 / s))
  res_fisher[, 1, 2] <- ((m * p^2 - m * p) * (1 / (m * s + 1))^(1 / s) * s - (p^2 + (m * p^2 - m * p) * s - p) * (1 / (m * s + 1))^(1 / s)*log(m * s + 1))/(m^2 * p * s^4 + 2 * m * p * s^3 + p * s^2 - ((m^2 * p - m^2) * s^4 + 2 * (m * p - m) * s^3 + (p - 1) * s^2) * (1 / (m * s + 1))^(1 / s))
  res_fisher[, 1, 3] <- -(1 / (m * s + 1))^(1 / s)/(m * p * s - ((m * p - m) * s + p - 1) * (1 / (m * s + 1))^(1 / s) + p)
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (log(1 + s * m) + digamma(1 / s) - digamma(m + 1 / s))^2 / s^4 * (1 - p - (1 - p) * (1 + s * m)^(-1 / s)) + ((1 - p)^2 * ((1 + s * m) * log(1 + s * m) - s * m)) / (s^4 * (1 + s * m)^(2 + 1 / s) * (1 + p * (1 + s * m)^(1 / s) - p))
  res_fisher[, 2, 3] <- -((1 / (m * s + 1))^(1 / s) * m * s - (m * s + 1) * (1 / (m * s + 1))^(1 / s) * log(m * s + 1)) / (m * p * s^3 + p * s^2 - ((m * p - m) * s^3 + (p - 1) * s^2) * (1 / (m * s + 1))^(1 / s))
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 3, 2] <- res_fisher[, 3, 2]
  res_fisher[, 3, 3] <- ((1 / (m * s + 1))^(1 / s) - 1) / (p^2 - (p^2 - 2 * p + 1) * (1 / (m * s + 1))^(1 / s) - p)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_zinegbin_nb2_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  p <- f[3]
  res_random <- sample(c(0L, NA_real_), size = t, replace = TRUE, prob = c(p, 1 - p))
  res_random[is.na(res_random)] <- be_silent(stats::rnbinom(sum(is.na(res_random)), mu = m, size = 1 / s))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_zinegbin_nb2_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_zero <- mean(y == 0L, na.rm = TRUE)
  p <- 0
  m <- y_mean
  s <- (y_var - y_mean) / y_mean^2
  for (i in 1:1e3) {
    p <- (y_zero - (1 / (1 + s * m))^(1 / s)) / (1 - (1 / (1 + s * m))^(1 / s))
    m <- y_mean / (1 - p)
    s <- y_var / y_mean / m - 1 / m - p
  }
  p <- max(min(p, 1 - 1e-6), 1e-6)
  m <- max(m, 1e-6)
  s <- max(s, 1e-6)
  res_start <- c(m, s, p)
  return(res_start)
}
# ------------------------------------------------------------------------------


