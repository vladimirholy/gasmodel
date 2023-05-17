
# VON MISES DISTRIBUTION / MEAN-CONCENTRATION PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_vonmises_meanconc_parameters <- function(n) {
  group_of_par_names <- c("mean", "conc")
  par_names <- c("mean", "conc")
  par_support <- c("real", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_vonmises_meanconc_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_density <- be_silent(CircStats::dvm(y, mu = m, kappa = v))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_vonmises_meanconc_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_loglik <- be_silent(log(CircStats::dvm(y, mu = m, kappa = v)))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_vonmises_meanconc_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_vonmises_meanconc_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_var <- 1 - besselI(v, nu = 1) / besselI(v, nu = 0)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_vonmises_meanconc_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- v * sin(y - m)
  res_score[, 2] <- cos(y - m) - besselI(v, nu = 1) / besselI(v, nu = 0)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_vonmises_meanconc_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- v * besselI(v, nu = 1) / besselI(v, nu = 0)
  res_fisher[, 1, 2] <- 0
  res_fisher[, 2, 1] <- 0
  res_fisher[, 2, 2] <- 1 / 2 - besselI(v, nu = 1)^2 / besselI(v, nu = 0)^2 + besselI(v, nu = 2) / (2 * besselI(v, nu = 0))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_vonmises_meanconc_random <- function(t, f) {
  m <- f[1]
  v <- f[2]
  res_random <- be_silent(CircStats::rvm(t, mean = m, k = v))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_vonmises_meanconc_start <- function(y) {
  ml_est <- CircStats::vm.ml(y)
  m <- ml_est[1, 1] + (2 * pi) * (ml_est[1, 1] < 0)
  v <- ml_est[1, 2]
  res_start <- c(m, v)
  return(res_start)
}
# ------------------------------------------------------------------------------


