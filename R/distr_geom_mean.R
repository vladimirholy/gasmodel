
# GEOMETRIC DISTRIBUTION / MEAN PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_geom_mean_parameters <- function(n) {
  group_of_par_names <- c("mean")
  par_names <- c("mean")
  par_support <- c("positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_geom_mean_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  res_density <- be_silent(stats::dgeom(y, prob = 1 / (1 + m)))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_geom_mean_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  res_loglik <- be_silent(stats::dgeom(y, prob = 1 / (1 + m), log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_geom_mean_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_geom_mean_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  res_var <- m * (1 + m)
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_geom_mean_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 1L)
  res_score[, 1] <- (y - m) / (m^2 + m)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_geom_mean_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 1L, 1L))
  res_fisher[, 1, 1] <- 1 / (m^2 + m)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_geom_mean_random <- function(t, f) {
  m <- f[1]
  res_random <- be_silent(stats::rgeom(t, prob = 1 / (1 + m)))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_geom_mean_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  m <- max(y_mean, 1e-6)
  res_start <- m
  return(res_start)
}
# ------------------------------------------------------------------------------


