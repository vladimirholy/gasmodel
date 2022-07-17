
# MULTIVARIATE STUDENTS T DISTRIBUTION / MEAN-VARIANCE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_mvt_meanvar_parameters <- function(n) {
  group_of_par_names <- c(rep("mean", n), rep("var", n), rep("cov", n * (n - 1) / 2), "df")
  par_names <- c(paste0("mean", 1:n), paste0("var", 1:n), paste0("cov", unlist(lapply(1:(n - 1), function(i) { return(paste0(i, ".", (i + 1):n)) }))), "df")
  par_support <- c(rep("real", n), rep("positive", n), rep("real", n * (n - 1) / 2), "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_mvt_meanvar_density <- function(y, f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * (ncol(f) - 1)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  v <- f[, 2 * n + n * (n - 1) / 2 + 1]
  res_density <- matrix(0, nrow = t, ncol = 1L)
  for (i in 1:t) {
    res_density[i, ] <- suppressWarnings(mvnfast::dmvt(y[i, ], mu = m[i, ], sigma = convert_varcov_vector_to_varcov_matrix(sc[i, ]), df = v[i], log = FALSE))
  }
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_mvt_meanvar_loglik <- function(y, f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * (ncol(f) - 1)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  v <- f[, 2 * n + n * (n - 1) / 2 + 1]
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  for (i in 1:t) {
    res_loglik[i, ] <- suppressWarnings(mvnfast::dmvt(y[i, ], mu = m[i, ], sigma = convert_varcov_vector_to_varcov_matrix(sc[i, ]), df = v[i], log = TRUE))
  }
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_mvt_meanvar_mean <- function(f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * (ncol(f) - 1)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  v <- f[, 2 * n + n * (n - 1) / 2 + 1]
  res_mean <- m
  res_mean[v <= 1, ] <- NA_real_
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_mvt_meanvar_var <- function(f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * (ncol(f) - 1)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  v <- f[, 2 * n + n * (n - 1) / 2 + 1]
  res_var <- array(0, dim = c(t, n, n))
  for (i in 1:t) {
    res_var[i, , ] <- convert_varcov_vector_to_varcov_matrix(sc[i, ]) * v[i] / (v[i] - 2)
  }
  res_var[v <= 2, , ] <- NA_real_
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_mvt_meanvar_score <- function(y, f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * (ncol(f) - 1)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  v <- f[, 2 * n + n * (n - 1) / 2 + 1]
  res_score <- matrix(0, nrow = t, ncol = 2 * n + n * (n - 1) / 2 + 1)
  for (i in 1:t) {
    sigma_inv <- matrix_inv(convert_varcov_vector_to_varcov_matrix(sc[i, ]))
    sigma_ysy <- as.vector((y[i, ] - m[i, ]) %*% sigma_inv %*% (y[i, ] - m[i, ]))
    psi_inv <- convert_varcov_matrix_to_sc_vector(sigma_inv)
    psi_syys <- convert_varcov_matrix_to_sc_vector(sigma_inv %*% (y[i, ] - m[i, ]) %*% (y[i, ] - m[i, ]) %*% sigma_inv)
    res_score[i, 1:n] <- (v[i] + n) / (v[i] + sigma_ysy) * sigma_inv %*% (y[i, ] - m[i, ])
    res_score[i, (n + 1):(2 * n + n * (n - 1) / 2)] <- 1 / 2 * (v[i] + n) / (v[i] + sigma_ysy) * psi_syys - 1 / 2 * psi_inv
    res_score[i, 2 * n + n * (n - 1) / 2 + 1] <- 1 / 2 * (sigma_ysy  - n) / (sigma_ysy + v[i]) - 1 / 2 * log(1 + 1 / v[i] * sigma_ysy) - 1 / 2 * digamma(v[i] / 2) + 1 / 2 * digamma((v[i] + n) / 2)
  }
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_mvt_meanvar_fisher <- function(f) {
  t <- nrow(f)
  n <- sqrt(9 / 4 + 2 * (ncol(f) - 1)) - 3 / 2
  m <- f[, 1:n, drop = FALSE]
  sc <- f[, (n + 1):(2 * n + n * (n - 1) / 2), drop = FALSE]
  v <- f[, 2 * n + n * (n - 1) / 2 + 1]
  res_fisher <- array(0, dim = c(t, 2 * n + n * (n - 1) / 2 + 1, 2 * n + n * (n - 1) / 2 + 1))
  for (i in 1:t) {
    sigma_inv <- matrix_inv(convert_varcov_vector_to_varcov_matrix(sc[i, ]))
    psi_inv <- convert_varcov_matrix_to_sc_vector(sigma_inv)
    psi_krone <- convert_krone_matrix_to_sc_matrix(kronecker(sigma_inv, sigma_inv))
    res_fisher[i, 1:n, 1:n] <- (v[i] + n) / (v[i] + n + 2) * sigma_inv
    res_fisher[i, (n + 1):(2 * n + n * (n - 1) / 2), (n + 1):(2 * n + n * (n - 1) / 2)] <- 1 / 4 * (v[i] + n) / (v[i] + n + 2) * psi_krone + 1 / 4 * (v[i] + n - 2) / (v[i] + n + 2) * psi_inv %*% t(psi_inv)
    res_fisher[i, (n + 1):(2 * n + n * (n - 1) / 2), 2 * n + n * (n - 1) / 2 + 1] <-  -1 / ((v[i] + n + 2) * (v[i] + n)) * psi_inv
    res_fisher[i, 2 * n + n * (n - 1) / 2 + 1, (n + 1):(2 * n + n * (n - 1) / 2)] <- res_fisher[i, (n + 1):(2 * n + n * (n - 1) / 2), 2 * n + n * (n - 1) / 2 + 1]
    res_fisher[i, 2 * n + n * (n - 1) / 2 + 1, 2 * n + n * (n - 1) / 2 + 1] <-  -1 / 2 * n * (n + v[i] + 4) / (v[i] * (n + v[i]) * (n + v[i] + 2)) + 1 / 4 * trigamma(v[i] / 2) - 1 / 4 * trigamma((v[i] + n) / 2)
  }
  return(res_fisher)
}
# ------------------------------------------------------------------------------

# Random Generation Function ---------------------------------------------------

distr_mvt_meanvar_random <- function(t, f) {
  n <- sqrt(9 / 4 + 2 * (length(f) - 1)) - 3 / 2
  m <- f[1:n]
  sc <- f[(n + 1):(2 * n + n * (n - 1) / 2)]
  v <- f[2 * n + n * (n - 1) / 2 + 1]
  res_random <- suppressWarnings(mvnfast::rmvt(t, mu = m, sigma = convert_varcov_vector_to_varcov_matrix(sc), df = v))
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_mvt_meanvar_start <- function(y) {
  y_mean <- colMeans(y, na.rm = TRUE)
  y_var <- stats::cov(y, use = "na.or.complete")
  y_kurt <- apply(y, MARGIN = 2, FUN = function(yy) { mean((yy - mean(yy, na.rm = TRUE))^4) / stats::var(yy, na.rm = TRUE)^2 - 3 })
  m <- y_mean
  sc <- convert_varcov_matrix_to_varcov_vector(y_var * (3 + mean(y_kurt)) / (3 + 2 * mean(y_kurt)))
  sc[1:ncol(y)] <- pmax(sc[1:ncol(y)], 1e-6)
  message(y_kurt)
  v <- max(6 / mean(y_kurt) + 4, 2.1)
  res_start <- c(m, sc, v)
  return(res_start)
}
# ------------------------------------------------------------------------------


