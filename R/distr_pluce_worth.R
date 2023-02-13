
# PLACKETT-LUCE DISTRIBUTION / WORTH PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_pluce_worth_parameters <- function(n) {
  group_of_par_names <- rep("worth", times = n)
  par_names <- paste0("worth", 1:n)
  par_support <- rep("positive", times = n)
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_pluce_worth_density <- function(y, f) {
  res_density <- exp(distr_pluce_worth_loglik(y = y, f = f))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_pluce_worth_loglik <- function(y, f) {
  t <- nrow(f)
  n <- ncol(f)
  res_loglik <- matrix(0, nrow = t, ncol = 1L)
  y_inv <- matrix(NA_real_, nrow = t, ncol = n)
  for (i in 1:t) {
    k <- sum(is.finite(y[i, ]))
    if (k < n) {
      y[i, !is.finite(y[i, ])] <- (k + 1):n
    }
    y_inv[i, ] <- Matrix::invPerm(y[i, ])
    o1 <- t(matrix(f[i, y_inv[i, ]], nrow = n, ncol = n))
    o1[lower.tri(o1)] <- 0
    o1 <- o1[1:k, 1:n]
    o2 <- log(rowSums(o1))
    res_loglik[i, ] <- sum(log(f[i, y_inv[i, 1:k]])) - sum(o2)
  }
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_pluce_worth_mean <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_mean <- matrix(0, nrow = t, ncol = n)
  if (n <= 6L) {
    y_all <- arrangements::permutations(n)
    for (i in 1:t) {
      for (j in 1:nrow(y_all)) {
        prob <- distr_pluce_worth_density(y = y_all[j, , drop = FALSE], f = f[i, , drop = FALSE])
        res_mean[i, ] <- res_mean[i, ] + y_all[j, ] * as.vector(prob)
      }
    }
  } else {
    if(!exists(".Random.seed")) { set.seed(NULL) }
    old_seed <- .Random.seed
    set.seed(13)
    for (i in 1:t) {
      if (all(is.finite(f[i, ])) && all(f[i, ] > 1e-12) && all(f[i, ] < 1e12)) {
        y_rand <- distr_pluce_worth_random(t = 1e3, f = f[i, , drop = FALSE])
        res_mean[i, ] <- colMeans(y_rand)
      } else {
        res_mean[i, ] <- NA_real_
      }
    }
    .Random.seed <- old_seed
  }
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_pluce_worth_var <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_var <- array(0, dim = c(t, n, n))
  if (n <= 6L) {
    supp_mean <- matrix(0, nrow = t, ncol = n)
    supp_square <- array(0, dim = c(t, n, n))
    y_all <- arrangements::permutations(n)
    for (i in 1:t) {
      for (j in 1:nrow(y_all)) {
        prob <- distr_pluce_worth_density(y = y_all[j, , drop = FALSE], f = f[i, , drop = FALSE])
        supp_mean[i, ] <- supp_mean[i, ] + y_all[j, ] * as.vector(prob)
        supp_square[i, , ] <- supp_square[i, , ] + t(y_all[j, , drop = FALSE]) %*% y_all[j, , drop = FALSE] * as.vector(prob)
      }
      res_var[i, , ] <- supp_square[i, , ] - t(supp_mean[i, , drop = FALSE]) %*% supp_mean[i, , drop = FALSE]
    }
  } else {
    if(!exists(".Random.seed")) { set.seed(NULL) }
    old_seed <- .Random.seed
    set.seed(13)
    for (i in 1:t) {
      if (all(is.finite(f[i, ])) && all(f[i, ] > 1e-12) && all(f[i, ] < 1e12)) {
        y_rand <- distr_pluce_worth_random(t = 1e3, f = f[i, , drop = FALSE])
        res_var[i, , ] <- stats::cov(y_rand)
      } else {
        res_var[i, , ] <- NA_real_
      }
    }
    .Random.seed <- old_seed
  }
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_pluce_worth_score <- function(y, f) {
  t <- nrow(f)
  n <- ncol(f)
  res_score <- matrix(0, nrow = t, ncol = n)
  y_inv <- matrix(NA_real_, nrow = t, ncol = n)
  for (i in 1:t) {
    k <- sum(is.finite(y[i, ]))
    if (k < n) {
      y[i, !is.finite(y[i, ])] <- (k + 1):n
    }
    y_inv[i, ] <- Matrix::invPerm(y[i, ])
    o1 <- t(matrix(f[i, y_inv[i, ]], nrow = n, ncol = n))
    o1[lower.tri(o1)] <- 0
    o1 <- o1[1:k, 1:n]
    o2 <- 1 / rowSums(o1)
    o3 <- t(matrix(o2, nrow = k, ncol = n))
    o3[upper.tri(o3)] <- 0
    res_score[i, ] <- (y[i, ] <= k) / f[i, ] - rowSums(o3[y[i, ], ])
  }
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_pluce_worth_fisher <- function(f) {
  t <- nrow(f)
  n <- ncol(f)
  res_fisher <- array(0, dim = c(t, n, n))
  if (n <= 6L) {
    y_all <- arrangements::permutations(n)
    for (i in 1:t) {
      for (j in 1:nrow(y_all)) {
        prob <- distr_pluce_worth_density(y = y_all[j, , drop = FALSE], f = f[i, , drop = FALSE])
        score <- distr_pluce_worth_score(y = y_all[j, , drop = FALSE], f = f[i, , drop = FALSE])
        res_fisher[i, , ] <- res_fisher[i, , ] + t(score) %*% score * as.vector(prob)
      }
    }
  } else {
    if(!exists(".Random.seed")) { set.seed(NULL) }
    old_seed <- .Random.seed
    set.seed(13)
    for (i in 1:t) {
      if (all(is.finite(f[i, ])) && all(f[i, ] > 1e-12) && all(f[i, ] < 1e12)) {
        y_rand <- distr_pluce_worth_random(t = 1e3, f = f[i, , drop = FALSE])
        for (j in 1:nrow(y_rand)) {
          score <- distr_pluce_worth_score(y = y_rand[j, , drop = FALSE], f = f[i, , drop = FALSE])
          res_fisher[i, , ] <- res_fisher[i, , ] + t(score) %*% score / nrow(y_rand)
        }
      } else {
        res_fisher[i, , ] <- NA_real_
      }
    }
    .Random.seed <- old_seed
  }
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_pluce_worth_random <- function(t, f) {
  n <- length(f)
  f <- pmax(pmin(1e12, f), 1e-12)
  f[is.na(f)] <- 1
  w <- f / prod(f)^(1 / n)
  res_random <- t(replicate(Matrix::invPerm(sample(1:n, replace = FALSE, prob = w)), n = t))
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_pluce_worth_start <- function(y) {
  y <- y[stats::complete.cases(y), ]
  t <- nrow(y)
  n <- ncol(y)
  my_fn <- function(f_tilde) { mean(as.vector(distr_pluce_worth_loglik(y = y, f = matrix(exp(c(f_tilde, 0)), nrow = t, ncol = n, byrow = TRUE)))) }
  my_gr <- function(f_tilde) { colMeans(distr_pluce_worth_score(y = y, f = matrix(exp(c(f_tilde, 0)), nrow = t, ncol = n, byrow = TRUE)))[1:(n - 1)] * exp(f_tilde) }
  my_par <- rep(0, n - 1)
  my_optim <- stats::optim(par = my_par, fn = my_fn, gr = my_gr, method = "BFGS", control = list(fnscale = -1, maxit = 100))
  res_start <- exp(c(my_optim$par, 0) - mean(c(my_optim$par, 0)))
  return(res_start)
}
# ------------------------------------------------------------------------------


