
# ZERO-INFLATED SKELLAM DISTRIBUTION / DIFFERENCE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_ziskellam_diff_parameters <- function(n) {
  group_of_par_names <- c("rate", "rate", "infl")
  par_names <- c("rate1", "rate2", "infl")
  par_support <- c("positive", "positive", "probability")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_ziskellam_diff_density <- function(y, f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_density <- be_silent((y == 0L) * p + (1 - p) * exp(-(r1 + r2)) * (r1 / r2)^(y / 2) * besselI(x = 2 * sqrt(r1 * r2), nu = y))
  res_density[!is.finite(res_density)] <- -Inf
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_ziskellam_diff_loglik <- function(y, f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  res_loglik[y != 0L, ] <- be_silent(log(1 - p[y != 0L, ]) + y[y != 0L, ] / 2 * log(r1[y != 0L, ] / r2[y != 0L, ]) - r1[y != 0L, ] - r2[y != 0L, ] + log(besselI(x = 2 * sqrt(r1[y != 0L, ] * r2[y != 0L, ]), nu = y[y != 0L, ])))
  res_loglik[y == 0L, ] <- be_silent(log(p[y == 0L, ] + (1 - p[y == 0L, ]) * exp(-(r1[y == 0L, ] + r2[y == 0L, ])) * besselI(x = 2 * sqrt(r1[y == 0L, ] * r2[y == 0L, ]), nu = 0)))
  res_loglik[!is.finite(res_loglik)] <- -Inf
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_ziskellam_diff_mean <- function(f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_mean <- (1 - p) * (r1 - r2)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_ziskellam_diff_var <- function(f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  res_var <- (1 - p) * (p * (r1 - r2)^2 + r1 + r2)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_ziskellam_diff_score <- function(y, f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  bi_0 <- besselI(x = 2 * sqrt(r1 * r2), nu = 0)
  bi_1 <- besselI(x = 2 * sqrt(r1 * r2), nu = 1)
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- (y == 0L) * (((p - 1) * (sqrt(r1 * r2) * bi_0 - r2 * bi_1)) / (sqrt(r1 * r2) * (p * exp(r1 + r2) + (1 - p) * bi_0))) + (y != 0L) * (sqrt(r2 / r1) * besselI(x = 2 * sqrt(r1 * r2), nu = y - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = y) - 1)
  res_score[, 2] <- (y == 0L) * (((p - 1) * (sqrt(r1 * r2) * bi_0 - r1 * bi_1)) / (sqrt(r1 * r2) * (p * exp(r1 + r2) + (1 - p) * bi_0))) + (y != 0L) * (sqrt(r1 / r2) * besselI(x = 2 * sqrt(r1 * r2), nu = y - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = y) - y / r2 - 1)
  res_score[, 3] <- (y == 0L) * ((exp(r1 + r2) - bi_0) / (p * exp(r1 + r2) + (1 - p) * bi_0)) + (y != 0L) * (1 / (p - 1))
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_ziskellam_diff_fisher <- function(f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  p <- f[, 3, drop = FALSE]
  bi_0 <- besselI(x = 2 * sqrt(r1 * r2), nu = 0)
  bi_1 <- besselI(x = 2 * sqrt(r1 * r2), nu = 0)
  rat_bi_r <- besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2)
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  res_fisher[, 1, 1] <- (1 - p) * (1 - exp(-r1 - r2) * bi_0) * (1 - sqrt(r2 / r1) * rat_bi_r)^2 + ((1 - p)^2 * exp(-r1 - r2) * (sqrt(r1 * r2) * bi_0 - r2 * bi_1)^2) / (r1 * r2 * (p * exp(r1 + r2) + (1 - p) * bi_0))
  res_fisher[, 1, 2] <- (1 - p) * (1 - exp(-r1 - r2) * bi_0) * (1 - sqrt(r2 / r1) * rat_bi_r) * (r1 / r2 - sqrt(r1 / r2) * rat_bi_r) + ((1 - p)^2 * exp(-r1 - r2) * (sqrt(r1 * r2) * bi_0 - r2 * bi_1) * (sqrt(r1 * r2) * bi_0 - r1 * bi_1)) / (r1 * r2 * (p * exp(r1 + r2) + (1 - p) * bi_0))
  res_fisher[, 1, 3] <- ((p - 1) * (1 - exp(-r1 - r2) * bi_0) * (sqrt(r1 * r2) * bi_0 - r2 * bi_1)) / (sqrt(r1 * r2) * (p * exp(r1 + r2) + (1 - p) * bi_0))
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (1 - p) * (1 - exp(-r1 - r2) * bi_0) * (r1 / r2 - sqrt(r1 / r2) * rat_bi_r)^2 + ((1 - p)^2 * exp(-r1 - r2) * (sqrt(r1 * r2) * bi_0 - r1 * bi_1)^2) / (r1 * r2 * (p * exp(r1 + r2) + (1 - p) * bi_0))
  res_fisher[, 2, 3] <- ((p - 1) * (1 - exp(-r1 - r2) * bi_0) * (sqrt(r1 * r2) * bi_0 - r1 * bi_1)) / (sqrt(r1 * r2) * (p * exp(r1 + r2) + (1 - p) * bi_0))
  res_fisher[, 3, 1] <- res_fisher[, 1, 3]
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <- (exp(r1 + r2) - bi_0) / ((1 - p) * (p * exp(r1 + r2) + (1 - p) * bi_0))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_ziskellam_diff_random <- function(t, f) {
  r1 <- f[1]
  r2 <- f[2]
  p <- f[3]
  res_random <- sample(c(0L, NA_real_), size = t, replace = TRUE, prob = c(p, 1 - p))
  res_random[is.na(res_random)] <- be_silent(stats::rpois(sum(is.na(res_random)), lambda = r1) - stats::rpois(sum(is.na(res_random)), lambda = r2))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_ziskellam_diff_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_zero <- mean(y == 0L, na.rm = TRUE)
  p <- 0
  r1 <- (y_var + y_mean) / 2
  r2 <- (y_var - y_mean) / 2
  for (i in 1:1e3) {
    p <- (y_zero - exp(-(r1 + r2)) * besselI(x = 2 * sqrt(r1 * r2), nu = y)) / (1 - exp(-(r1 + r2)) * besselI(x = 2 * sqrt(r1 * r2), nu = y))
    r1 <- y_mean / (1 - p) + r2
    r2 <- y_var / (1 - p) - p * y_mean^2 / (1 - p)^2 - r1
  }
  p <- max(min(p, 1 - 1e-6), 1e-6)
  r1 <- max(r1, 1e-6)
  r2 <- max(r2, 1e-6)
  res_start <- c(r1, r2, p)
  return(res_start)
}
# ------------------------------------------------------------------------------


