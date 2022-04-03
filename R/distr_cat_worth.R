
# CATEGORICAL DISTRIBUTION / WORTH PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_cat_worth_parameters <- function(n) {
  group_of_par_names <- rep("worth", times = n)
  par_names <- paste0("worth", 1:n)
  par_support <- rep("positive", times = n)
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_cat_worth_density <- function(y, f, n) {
  t <- nrow(f)
  res_density <- f[y == 1L] / rowSums(f)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_cat_worth_loglik <- function(y, f, supp) {
  t <- nrow(f)
  res_loglik <- log(f[y == 1L]) - log(rowSums(f))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_cat_worth_mean <- function(f) {
  res_mean <- f / rowSums(f)
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_cat_worth_var <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_var <- array(0, dim = c(t, n, n))
  for (i in 1:t) {
    res_var[i, , ] <- diag(f[i, ] / sum(f[i, ])) - f[i, ] %*% t(f[i, ]) / sum(f[i, ])^2
  }
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_cat_worth_score <- function(y, f) {
  t <- nrow(f)
  n <- ncol(f)
  res_score <- matrix(0, nrow = t, ncol = n)
  for (i in 1:t) {
    res_score[i, ] <- y[i, ] / f[i, ] - 1 / sum(f[i, ])
  }
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_cat_worth_fisher <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_fisher <- array(0, dim = c(t, n, n))
  for (i in 1:t) {
    res_fisher[i, , ] <- diag(sum(f[i, ]) / f[i, ]) - 1 / sum(f[i, ])^2
  }
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_cat_worth_random <- function(t, f) {
  n <- length(f)
  res_random <- 1 * (matrix(1:n, nrow = t, ncol = n, byrow = TRUE) == sample(1:n, size = t, replace = TRUE, prob = f))
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_cat_worth_start <- function(y) {
  y <- y[stats::complete.cases(y), ]
  res_start <- colMeans(y)
  return(res_start)
}
# ------------------------------------------------------------------------------


