
# SKELLAM DISTRIBUTION / DIFFERENCE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_skellam_diff_parameters <- function(n) {
  group_of_par_names <- c("rate", "rate")
  par_names <- c("rate1", "rate2")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_skellam_diff_density <- function(y, f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  res_density <- be_silent(exp(-(r1 + r2)) * (r1 / r2)^(y / 2) * besselI(x = 2 * sqrt(r1 * r2), nu = y))
  res_density[!is.finite(res_density)] <- -Inf
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_skellam_diff_loglik <- function(y, f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  res_loglik <- be_silent(y / 2 * log(r1 / r2) - r1 - r2 + log(besselI(x = 2 * sqrt(r1 * r2), nu = y)))
  res_loglik[!is.finite(res_loglik)] <- -Inf
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_skellam_diff_mean <- function(f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  res_mean <- r1 - r2
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_skellam_diff_var <- function(f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  res_var <- r1 + r2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_skellam_diff_score <- function(y, f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- sqrt(r2 / r1) * besselI(x = 2 * sqrt(r1 * r2), nu = y - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = y) - 1
  res_score[, 2] <- sqrt(r1 / r2) * besselI(x = 2 * sqrt(r1 * r2), nu = y - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = y) - y / r2 - 1
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_skellam_diff_fisher <- function(f) {
  t <- nrow(f)
  r1 <- f[, 1, drop = FALSE]
  r2 <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- r2 / r1 * (besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2))^2 - 2 * sqrt(r2 / r1) * besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2) + 1
  res_fisher[, 1, 2] <- (besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2))^2 - 2 * sqrt(r1 / r2) * besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2) + r1 / r2
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- r1 / r2 * (besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2))^2 - 2 * (r1 / r2)^(3 / 2) * besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2 - 1) / besselI(x = 2 * sqrt(r1 * r2), nu = r1 - r2) + (r1 / r2)^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_skellam_diff_random <- function(t, f) {
  r1 <- f[1]
  r2 <- f[2]
  res_random <- be_silent(stats::rpois(t, lambda = r1) - stats::rpois(t, lambda = r2))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_skellam_diff_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  r1 <- max((y_var + y_mean) / 2, 1e-6)
  r2 <- max((y_var - y_mean) / 2, 1e-6)
  res_start <- c(r1, r2)
  return(res_start)
}
# ------------------------------------------------------------------------------


