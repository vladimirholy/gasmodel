
# ZERO-INFLATED SKELLAM DISTRIBUTION / MEAN-VARIANCE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_ziskellam_meanvar_parameters <- function(n) {
  group_of_par_names <- c("mean", "var", "infl")
  par_names <- c("mean", "var", "infl")
  par_support <- c("real", "positive", "probability")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_ziskellam_meanvar_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  p[s <= abs(m)] <- NA_real_
  res_density <- suppressWarnings((y == 0L) * p + (1 - p) * exp(-s) * ((s + m) / (s - m))^(y / 2) * besselI(x = sqrt(s^2 - m^2), nu = y))
  res_density[!is.finite(res_density)] <- -Inf
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_ziskellam_meanvar_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  p[s <= abs(m)] <- NA_real_
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  res_loglik[y != 0L, ] <- suppressWarnings(log(1 - p[y != 0L, ]) + y[y != 0L, ] / 2 * log((s[y != 0L, ] + m[y != 0L, ]) / (s[y != 0L, ] - m[y != 0L, ])) - s[y != 0L, ] + log(besselI(x = sqrt(s[y != 0L, ]^2 - m[y != 0L, ]^2), nu = y[y != 0L, ])))
  res_loglik[y == 0L, ] <- suppressWarnings(log(p[y == 0L, ] + (1 - p[y == 0L, ]) * exp(-s[y == 0L, ]) * besselI(x = sqrt(s[y == 0L, ]^2 - m[y == 0L, ]^2), nu = 0)))
  res_loglik[!is.finite(res_loglik)] <- -Inf
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_ziskellam_meanvar_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  p[s <= abs(m)] <- NA_real_
  res_mean <- (1 - p) * m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_ziskellam_meanvar_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  p[s <= abs(m)] <- NA_real_
  res_var <- (1 - p) * (s + p * m^2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_ziskellam_meanvar_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  p[s <= abs(m)] <- NA_real_
  bi_0 <- besselI(x = sqrt(s^2 - m^2), nu = 0)
  bi_1 <- besselI(x = sqrt(s^2 - m^2), nu = 1)
  tri_bi_y <- (besselI(x = sqrt(s^2 - m^2), nu = y - 1) + besselI(x = sqrt(s^2 - m^2), nu = y + 1)) / besselI(x = sqrt(s^2 - m^2), nu = y)
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- (y == 0L) * ((m * (p - 1) * bi_1) / (sqrt(s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))) + (y != 0L) * ((s * y) / (s^2 - m^2) - m / (2 * sqrt(s^2 - m^2)) * tri_bi_y)
  res_score[, 2] <- (y == 0L) * (((p - 1) * (sqrt(s^2 - m^2) * bi_0 - s * bi_1)) / (sqrt(s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))) + (y != 0L) * (-(m * y) / (s^2 - m^2) + s / (2 * sqrt(s^2 - m^2)) * tri_bi_y - 1)
  res_score[, 3] <- (y == 0L) * ((exp(s) - bi_0) / (p * exp(s) + (1 - p) * bi_0)) + (y != 0L) * (1 / (p - 1))
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_ziskellam_meanvar_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  p[s <= abs(m)] <- NA_real_
  bi_0 <- besselI(x = sqrt(s^2 - m^2), nu = 0)
  bi_1 <- besselI(x = sqrt(s^2 - m^2), nu = 1)
  tri_bi_m <- (besselI(x = sqrt(s^2 - m^2), nu = m - 1) + besselI(x = sqrt(s^2 - m^2), nu = m + 1)) / besselI(x = sqrt(s^2 - m^2), nu = m)
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- (m^2 * (1 - p) * (1 - exp(-s) * bi_0)) / (4 * (s^2 - m^2)) * ((2 * s) / (sqrt(s^2 - m^2)) - tri_bi_m)^2 + (m^2 * (1 - p)^2 * exp(-s) * bi_1^2) / ((s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))
  res_fisher[, 1, 2] <- (m * s * (p - 1) * (1 - exp(-s) * bi_0)) / (4 * (s^2 - m^2)) * ((2 * s) / (sqrt(s^2 - m^2)) - tri_bi_m)^2 + (m * (1 - p)^2 * exp(-s) * bi_1 * (sqrt(s^2 - m^2) * bi_0 + s * bi_1)) / ((s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))
  res_fisher[, 1, 3] <- (m * (p - 1) * (1 - exp(-s) * bi_0) * bi_1) / (sqrt(s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (s^2 * (1 - p) * (1 - exp(-s) * bi_0)) / (4 * (s^2 - m^2)) * ((2 * s) / (sqrt(s^2 - m^2)) - tri_bi_m)^2 + ((1 - p)^2 * exp(-s) * (sqrt(s^2 - m^2) * bi_0 + s * bi_1)^2) / ((s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))
  res_fisher[, 2, 3] <- ((p - 1) * (1 - exp(-s) * bi_0) * (sqrt(s^2 - m^2) * bi_0 - s * bi_1)) / (sqrt(s^2 - m^2) * (p * exp(s) + (1 - p) * bi_0))
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <- (exp(s) - bi_0) / ((1 - p) * (p * exp(s) + (1 - p) * bi_0))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_ziskellam_meanvar_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  p <- f[3]
  res_random <- sample(c(0L, NA_real_), size = t, replace = TRUE, prob = c(p, 1 - p))
  res_random[is.na(res_random)] <- suppressWarnings(stats::rpois(sum(is.na(res_random)), lambda = (s + m) / 2) - stats::rpois(sum(is.na(res_random)), lambda = (s - m) / 2))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_ziskellam_meanvar_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_zero <- mean(y == 0L, na.rm = TRUE)
  p <- 0
  m <- y_mean
  s <- y_var
  for (i in 1:1e3) {
    p <- (y_zero - exp(-s) * besselI(x = sqrt(s^2 - m^2), nu = 0)) / (1 - exp(-s) * besselI(x = sqrt(s^2 - m^2), nu = 0))
    m <- y_mean / (1 - p)
    s <- y_var / (1 - p) - m^2 * p
  }
  p <- max(min(p, 1 - 1e-6), 1e-6)
  s <- max(y_var, abs(y_mean) + 1e-6)
  res_start <- c(m, s, p)
  return(res_start)
}
# ------------------------------------------------------------------------------


