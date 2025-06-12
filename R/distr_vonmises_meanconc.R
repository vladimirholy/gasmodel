
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
  res_density <- exp(v * cos(y - m)) / (2 * pi * besselI(v, nu = 0))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_vonmises_meanconc_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  v <- f[, 2, drop = FALSE]
  res_loglik <- v * cos(y - m) - log(2 * pi * besselI(v, nu = 0))
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
  res_fisher[, 2, 2] <- 1 / 2 - besselI(v, nu = 1)^2 / besselI(v, nu = 0)^2 + besselI(v, nu = 2) / (2 * besselI(v, nu = 0))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_vonmises_meanconc_random <- function(t, f) {
  m <- f[1]
  v <- f[2]
  res_random <- matrix(NA, nrow = t, ncol = 1L)
  a <- 1 + sqrt(1 + 4 * v^2)
  b <- (a - sqrt(2 * a)) / (2 * v)
  r <- (1 + b^2) / (2 * b)
  for (i in 1:t) {
    repeat {
      u1 <- stats::runif(1)
      z <- cos(pi * u1)
      f <- (1 + r * z) / (r + z)
      c <- v * (r - f)
      u2 <- stats::runif(1)
      if (u2 < c * (2 - c) || log(c / u2) + 1 - c >= 0) {
        u3 <- stats::runif(1)
        sign <- ifelse(u3 > 0.5, 1, -1)
        theta <- m + sign * acos(f)
        res_random[i, 1] <- theta %% (2 * pi)
        break()
      }
    }
  }
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_vonmises_meanconc_start <- function(y) {
  c <- mean(cos(y))
  s <- mean(sin(y))
  m <- atan2(s, c) %% (2 * pi)
  r <- sqrt(c^2 + s^2)
  if (r < 0.53) {
    v <- 2 * r + r^3 + (5 * r^5) / 6
  } else if (r < 0.85) {
    v <- -0.4 + 1.39 * r + 0.43 / (1 - r)
  } else {
    v <- 1 / (r^3 - 4 * r^2 + 3 * r)
  }
  res_start <- c(m, v)
  return(res_start)
}
# ------------------------------------------------------------------------------


