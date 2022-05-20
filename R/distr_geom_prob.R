
# GEOMETRIC DISTRIBUTION / PROBABILISTIC PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_geom_prob_parameters <- function(n) {
  group_of_par_names <- c("prob")
  par_names <- c("prob")
  par_support <- c("probability")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_geom_prob_density <- function(y, f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  res_density <- suppressWarnings(stats::dgeom(y, prob = p))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_geom_prob_loglik <- function(y, f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dgeom(y, prob = p, log = TRUE))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_geom_prob_mean <- function(f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  res_mean <- (1 - p) / p
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_geom_prob_var <- function(f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  res_var <- (1 - p) / p^2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_geom_prob_score <- function(y, f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 1L)
  res_score[, 1] <- (p * y + p - 1) / (p^2 - p)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_geom_prob_fisher <- function(f) {
  t <- nrow(f)
  p <- f[, 1, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 1L, 1L))
  res_fisher[, 1, 1] <- 1 / (p^2 - p^3)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_geom_prob_random <- function(t, f) {
  p <- f[1]
  res_random <- suppressWarnings(stats::rgeom(t, prob = p))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_geom_prob_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  p <- max(min(1 / (1 + y_mean), 1 - 1e-6), 1e-6)
  res_start <- p
  return(res_start)
}
# ------------------------------------------------------------------------------


