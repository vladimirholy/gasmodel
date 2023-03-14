
# SKELLAM DISTRIBUTION / MEAN-DISPERSION PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_skellam_meandisp_parameters <- function(n) {
  group_of_par_names <- c("mean", "disp")
  par_names <- c("mean", "disp")
  par_support <- c("real", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_skellam_meandisp_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_density <- be_silent(exp(-abs(m) - s) * ((abs(m) + m + s) / (abs(m) - m + s))^(y / 2) * besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = y))
  res_density[!is.finite(res_density)] <- -Inf
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_skellam_meandisp_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_loglik <- be_silent(y / 2 * log((abs(m) + m + s) / (abs(m) - m + s)) - abs(m) - s + log(besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = y)))
  res_loglik[!is.finite(res_loglik)] <- -Inf
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_skellam_meandisp_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_skellam_meandisp_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  res_var <- abs(m) + s
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_skellam_meandisp_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  tri_bi_y <- (besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = y - 1) + besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = y + 1)) / besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = y)
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- y / (2 * abs(m) + s) + (sign(m) * s) / (2 * sqrt(s^2 + 2 * abs(m) * s)) * tri_bi_y - sign(m)
  res_score[, 2] <- -(m * y) / (s^2 + 2 * abs(m) * s) + (abs(m) + s) / (2 * sqrt(s^2 + 2 * abs(m) * s)) * tri_bi_y - 1
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_skellam_meandisp_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  tri_bi_m <- (besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = m - 1) + besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = m + 1)) / besselI(x = sqrt(s^2 + 2 * abs(m) * s), nu = m)
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- s^2 / (4 * s^2 + 8 * abs(m) * s) * ((2 * abs(m) + 2 * s) / sqrt(s^2 + 2 * abs(m) * s) - tri_bi_m)^2
  res_fisher[, 1, 2] <- (sign(m) * (abs(m) + s) * s) / (4 * s^2 + 8 * abs(m) * s) * ((2 * abs(m) + 2 * s) / sqrt(s^2 + 2 * abs(m) * s) - tri_bi_m)^2
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- (abs(m) + s)^2 / (4 * s^2 + 8 * abs(m) * s) * ((2 * abs(m) + 2 * s) / sqrt(s^2 + 2 * abs(m) * s) - tri_bi_m)^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_skellam_meandisp_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  res_random <- be_silent(stats::rpois(t, lambda = (abs(m) + m + s) / 2) - stats::rpois(t, lambda = (abs(m) - m + s) / 2))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_skellam_meandisp_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  m <- y_mean
  s <- max(y_var - abs(y_mean), 1e-6)
  res_start <- c(m, s)
  return(res_start)
}
# ------------------------------------------------------------------------------


