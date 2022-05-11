
# STUDENTS T DISTRIBUTION / STANDARD PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_t_std_parameters <- function(n) {
  group_of_par_names <- c("location", "scale", "df")
  par_names <- c("location", "scale", "df")
  par_support <- c("real", "positive", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_t_std_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  v <- f[, 3, drop = FALSE]
  res_density <- suppressWarnings(stats::dt((y - m) / sqrt(s), df = v) / sqrt(s))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_t_std_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  v <- f[, 3, drop = FALSE]
  res_loglik <- suppressWarnings(stats::dt((y - m) / sqrt(s), df = v, log = TRUE) - log(s) / 2)
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_t_std_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  v <- f[, 3, drop = FALSE]
  res_mean <- m
  res_mean[v <= 1] <- Inf
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_t_std_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  v <- f[, 3, drop = FALSE]
  res_var <- s * v / (v - 2)
  res_var[v <= 2] <- Inf
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_t_std_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  v <- f[, 3, drop = FALSE]
  res_score <- matrix(0, nrow = t, ncol = 3L)
  res_score[, 1] <- (v + 1) * (y - m) / (m^2 - 2 * m * y + v * s + y^2)
  res_score[, 2] <- v * (m^2 - 2 * m * y - s + y^2) / (2 * s * (m^2 - 2 * m * y + v * s + y^2))
  res_score[, 3] <- m^2 / (2 * v^2 * s * ((y - m)^2 / (v * s) + 1)) + m^2 / (2 * v * s * ((y - m)^2 / (v * s) + 1)) + y^2 / (2 * v^2 * s * ((y - m)^2 / (v * s) + 1)) - (m * y) / (v^2 * s * ((y - m)^2 / (v * s) + 1)) + y^2 / (2 * v * s * ((y - m)^2 / (v * s) + 1)) - (m * y) / (v * s * ((y - m)^ 2 /(v * s) + 1)) - 1 / 2 * log((y - m)^2 / (v * s) + 1) - 1 / (2 * v) - digamma(v / 2) / 2 + digamma((v + 1) / 2) / 2
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_t_std_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  v <- f[, 3, drop = FALSE]
  res_fisher[, 1, 1] <- (v + 1) / ((v + 3) * s)
  res_fisher[, 2, 2] <- v / (2 * (v + 3) * s^2)
  res_fisher[, 2, 3] <- -1 / ((v + 1) * (v + 3) * s)
  res_fisher[, 3, 2] <- res_fisher[, 2, 3]
  res_fisher[, 3, 3] <-  1 / 4 * trigamma(v / 2) - 1 / 4 * trigamma((v + 1) / 2) - (v + 5) / (2 * v * (v + 1) * (v + 3))
  res_fisher <- array(0, dim = c(t, 3L, 3L))
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_t_std_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  v <- f[3]
  res_random <- suppressWarnings(m + sqrt(s) * stats::rt(t, df = v))
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_t_std_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  y_kurt <- mean((y - y_mean)^4) / y_var^2 - 3
  m <- y_mean
  s <- max(y_var * (3 + y_kurt) / (3 + 2 * y_kurt), 1e-6)
  v <- max(6 / y_kurt + 4, 2.1)
  res_start <- c(m, s, v)
  return(res_start)
}
# ------------------------------------------------------------------------------


