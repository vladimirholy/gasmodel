
# KUMARASWAMY DISTRIBUTION / CONCENTRATION PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_kuma_conc_parameters <- function(n) {
  group_of_par_names <- c("conc", "conc")
  par_names <- c("conc1", "conc2")
  par_support <- c("positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_kuma_conc_density <- function(y, f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_density <- a * b * y^(a - 1) * (1 - y^a)^(b - 1)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_kuma_conc_loglik <- function(y, f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_loglik <- log(a) + log(b) + (a - 1) * log(y) + (b - 1) * log(1 - y^a)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_kuma_conc_mean <- function(f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_mean <- b * beta(1 + 1 / a, b)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_kuma_conc_var <- function(f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_var <- b * beta(1 + 2 / a, b) - b^2 * beta(1 + 1 / a, b)^2
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_kuma_conc_score <- function(y, f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- ((b - 1) * log(y)) / (y^a - 1) + (a * b * log(y) + 1) / a
  res_score[, 2] <- log(1 - y^a) + 1 / b
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_kuma_conc_fisher <- function(f) {
  t <- nrow(f)
  a <- f[, 1, drop = FALSE]
  b <- f[, 2, drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- (1 + b / (b - 2) * ((digamma(b) - digamma(2))^2 - trigamma(b) + trigamma(2))) / a^2
  res_fisher[, 1, 2] <- -(digamma(b + 1) - digamma(2)) / (b - 1) / a
  res_fisher[, 2, 1] <- res_fisher[, 1, 2]
  res_fisher[, 2, 2] <- 1 / b^2
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_kuma_conc_random <- function(t, f) {
  a <- f[1]
  b <- f[2]
  res_random <- be_silent((1 - stats::runif(t)^(1 / b))^(1 / a))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_kuma_conc_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_square <- mean(y^2, na.rm = TRUE)
  a <- 1
  b <- 1
  for (i in 1:1e3) {
    b <- b + (y_square - b * beta(1 + 2 / a, b)) / ((b * digamma(b) - b * digamma(b + 2 / a + 1) + 1) * beta(1 + 2 / a, b))
    a <- a + (y_mean - b * beta(1 + 1 / a, b)) / (b * (digamma(b + 1 / a + 1) - digamma(1 + 1 / a)) * beta(1 + 1 / a, b) / a^2)
  }
  a <- max(a, 1e-6)
  b <- max(b, 1e-6)
  res_start <- c(a, b)
  return(res_start)
}
# ------------------------------------------------------------------------------


