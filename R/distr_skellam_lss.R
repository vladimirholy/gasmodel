
# POISSON DISTRIBUTION / LOCATION-SQ-SCALE PARAMETRIZATION


# Parameters Function ----------------------------------------------------------
distr_skellam_lss_parameters <- function(n) {
  group_of_par_names <- c("location", "sq.scale")
  par_names <- c("location", "sq.scale")
  par_support <- c("real", "positive")
  res_parameters <- list(group_of_par_names = group_of_par_names, par_names = par_names, par_support = par_support)
  return(res_parameters)
}
# ------------------------------------------------------------------------------


# Density Function -------------------------------------------------------------
distr_skellam_lss_density <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  res_density <- suppressWarnings(exp(-s) * ((s + m) / (s - m))^(y / 2) * besselI(x = sqrt(s^2 - m^2), nu = y))
  return(res_density)
}
# ------------------------------------------------------------------------------


# Log-Likelihood Function ------------------------------------------------------
distr_skellam_lss_loglik <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  res_loglik <- suppressWarnings(y / 2 * log((s + m) / (s - m)) - s + log(besselI(x = sqrt(s^2 - m^2), nu = y)))
  return(res_loglik)
}
# ------------------------------------------------------------------------------


# Mean Function ----------------------------------------------------------------
distr_skellam_lss_mean <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  res_mean <- m
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Variance Function ------------------------------------------------------------
distr_skellam_lss_var <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  res_var <- s
  res_var <- array(res_var, dim = c(t, 1, 1))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Score Function ---------------------------------------------------------------
distr_skellam_lss_score <- function(y, f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  res_score <- matrix(0, nrow = t, ncol = 2L)
  res_score[, 1] <- (s * y) / (s^2 - m^2) - m / (2 * sqrt(s^2 - m^2)) * (besselI(x = sqrt(s^2 - m^2), nu = y - 1) + besselI(x = sqrt(s^2 - m^2), nu = y + 1)) / besselI(x = sqrt(s^2 - m^2), nu = y)
  res_score[, 2] <- -(m * y) / (s^2 - m^2) + s / (2 * sqrt(s^2 - m^2)) * (besselI(x = sqrt(s^2 - m^2), nu = y - 1) + besselI(x = sqrt(s^2 - m^2), nu = y + 1)) / besselI(x = sqrt(s^2 - m^2), nu = y) - 1
  return(res_score)
}
# ------------------------------------------------------------------------------


# Fisher Information Function --------------------------------------------------
distr_skellam_lss_fisher <- function(f) {
  t <- nrow(f)
  m <- f[, 1, drop = FALSE]
  s <- f[, 2, drop = FALSE]
  m[s <= abs(m)] <- NA_real_
  s[s <= abs(m)] <- NA_real_
  res_fisher <- array(0, dim = c(t, 2L, 2L))
  res_fisher[, 1, 1] <- NA
  res_fisher[, 1, 2] <- NA
  res_fisher[, 2, 1] <- NA
  res_fisher[, 2, 2] <- NA
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Random Generation Function ---------------------------------------------------
distr_skellam_lss_random <- function(t, f) {
  m <- f[1]
  s <- f[2]
  if (s > abs(m)) {
    res_random <- suppressWarnings(stats::rpois(t, lambda = (v + m) / 2) - stats::rpois(t, lambda = (v - m) / 2))
  } else {
    res_random <- rep(NA_real_, times = t)
  }
  res_random <- matrix(res_random, nrow = t, ncol = 1L)
  return(res_random)
}
# ------------------------------------------------------------------------------


# Starting Estimates Function --------------------------------------------------
distr_skellam_lss_start <- function(y) {
  y_mean <- mean(y, na.rm = TRUE)
  y_var <- stats::var(y, na.rm = TRUE)
  m <- y_mean
  s <- max(y_var, abs(y_mean) + 1e-6)
  res_start <- c(m, s)
  return(res_start)
}
# ------------------------------------------------------------------------------


