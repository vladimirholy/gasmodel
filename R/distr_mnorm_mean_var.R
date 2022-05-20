
# MULTIVARIATE NORMAL DISTRIBUTION / MEAN-VARIANCE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_mnorm_mean_var_parameters <- function(n) {
  group_of_par_names <- c(rep("mean", n), rep("var", n), rep("cor", n * (n - 1) / 2))
  par_names <- c(paste0("mean", 1:n), paste0("var", 1:n), paste0("cor", unlist(lapply(1:(n - 1), function(i) { return(paste0(i, ".", (i + 1):n)) }))))
  par_support <- c(rep("real", n), rep("positive", n), rep("probability", n * (n - 1) / 2))
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_mnorm_mean_var_density <- function(y, f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * ncol(f)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  res_density <- matrix(0, nrow = t, ncol = 1L)
  for (i in 1:t) {
    res_density[i, ] <- suppressWarnings(mvnfast::dmvn(y[i, ], mu = m[i, ], sigma = convert_varcov_vector_to_varcov_matrix(sc[i, ])))
  }
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_mnorm_mean_var_loglik <- function(y, f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * ncol(f)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  for (i in 1:t) {
    res_loglik[i, ] <- suppressWarnings(mvnfast::dmvn(y[i, ], mu = m[i, ], sigma = convert_varcov_vector_to_varcov_matrix(sc[i, ]), log = TRUE))
  }
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_mnorm_mean_var_mean <- function(f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * ncol(f)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_mnorm_mean_var_var <- function(f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * ncol(f)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  res_var <- array(0, dim = c(t, n, n))
  for (i in 1:t) {
    res_var[i, , ] <- convert_varcov_vector_to_varcov_matrix(sc[i, ])
  }
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_mnorm_mean_var_score <- function(y, f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * ncol(f)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 2 * n + n * (n - 1) / 2)
  for (i in 1:t) {
    sigma_inv <- matrix_inv(convert_varcov_vector_to_varcov_matrix(sc[i, ]))
    sigma_der <- sigma_inv %*% (y[i, ] - m[i, ]) %*% (y[i, ] - m[i, ]) %*% sigma_inv - sigma_inv
    diag(sigma_der) <- diag(sigma_der) / 2
    res_score[i, 1:n] <- sigma_inv %*% (y[i, ] - m[i, ])
    res_score[i, (n + 1):(2 * n + n * (n - 1) / 2)] <- convert_varcov_matrix_to_varcov_vector(sigma_der)
  }
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_mnorm_mean_var_fisher <- function(f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * ncol(f)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  res_fisher <- array(0, dim = c(t, 2 * n + n * (n - 1) / 2, 2 * n + n * (n - 1) / 2))
  for (i in 1:t) {
    sigma_inv <- matrix_inv(convert_varcov_vector_to_varcov_matrix(sc[i, ]))
    sigma_fisher <- (kronecker(sigma_inv, sigma_inv) + matrix(as.vector(sigma_inv), ncol = 1) %*% matrix(as.vector(sigma_inv), nrow = 1)) / 4
    sigma_idx <- c(which(as.vector(as.logical(diag(n)))), which(as.vector(lower.tri(diag(n)))))
    res_fisher[i, 1:n, 1:n] <- sigma_inv
    res_fisher[i, (n + 1):(2 * n + n * (n - 1) / 2), (n + 1):(2 * n + n * (n - 1) / 2)] <- sigma_fisher[sigma_idx, sigma_idx]
  }
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_mnorm_mean_var_random <- function(t, f) {
  n <- sqrt(9 / 4 + 2 * length(f)) - 3 / 2
  m <- f[1:n]
  sc <- f[(n + 1):(2 * n + n * (n - 1) / 2)]
  res_random <- suppressWarnings(mvnfast::rmvn(t, mu = m, sigma = convert_varcov_vector_to_varcov_matrix(sc)))
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_mnorm_mean_var_start <- function(y) {
  y_mean <- colMeans(y, na.rm = TRUE)
  y_var <- stats::cov(y, use = "na.or.complete")
  m <- y_mean
  sc <- pmax(convert_varcov_matrix_to_varcov_vector(y_var), 1e-6)
  res_start <- c(m, sc)
  return(res_start)
}
# ------------------------------------------------------------------------------


