
# DIRICHLET DISTRIBUTION / CONCENTRATION PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_dirichlet_conc_parameters <- function(n) {
  group_of_par_names <- rep("conc", times = n)
  par_names <- paste0("conc", 1:n)
  par_support <- rep("positive", times = n)
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_dirichlet_conc_density <- function(y, f) {
  res_density <- exp(distr_dirichlet_conc_loglik(y = y, f = f))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_dirichlet_conc_loglik <- function(y, f) {
  t <- nrow(f)
  n <- ncol(f)
  res_loglik <- lgamma(rowSums(f)) + rowSums((f - 1) * log(y) - lgamma(f))
  res_loglik <- matrix(res_loglik, nrow = t, ncol = 1L)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_dirichlet_conc_mean <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_mean <- f / rowSums(f)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_dirichlet_conc_var <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  f_sum <- rowSums(f)
  res_mean <- f / f_sum
  res_var <- array(0, dim = c(t, n, n))
  for (i in 1:t) {
    res_var[i, , ] <- (diag(res_mean[i, ]) - t(res_mean[i, , drop = FALSE]) %*% res_mean[i, , drop = FALSE]) / (1 + f_sum[i])
  }
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_dirichlet_conc_score <- function(y, f) {
  t <- nrow(f)
  n <- ncol(f)
  res_score <- digamma(rowSums(f)) - digamma(f) + log(y)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_dirichlet_conc_fisher <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_fisher <- array(0, dim = c(t, n, n))
  for (i in 1:t) {
    res_fisher[i, , ] <- diag(trigamma(f[i, ])) - trigamma(sum(f[i, ]))
  }
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_dirichlet_conc_random <- function(t, f) {
  n <- length(f)
  res_random <- matrix(stats::rgamma(n * t, shape = f), nrow = t, ncol = n, byrow = TRUE)
  res_random <- res_random / rowSums(res_random)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_dirichlet_conc_start <- function(y) {
  y_mean <- colMeans(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  res_start <- pmax(y_mean * (y_mean * (1 - y_mean) / diag(y_var) - 1), 1e-6)
  return(res_start)
}
# ------------------------------------------------------------------------------


