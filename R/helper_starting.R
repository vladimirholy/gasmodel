
# STARTING ESTIMATE FUNCTION


# Estimate Starting Solution ---------------------------------------------------
starting_theta <- function(theta_start, theta_bound_lower, theta_bound_upper, data, model, fun, info_distr, info_par, info_coef, info_theta, print_progress) {
  par_start <- reparam_orig_to_tilde(f_orig = do.call(paste("distr", model$distr, model$param, "start", sep = "_"), args = list(y = data$y)), par_trans = info_par$par_trans)
  try_num <- 11L
  try_theta <- matrix(NA_real_, nrow = try_num, ncol = info_theta$theta_num)
  try_theta[1, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.10 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(0.01 / model$p[i], times = model$p[i]), phi = rep(0.10 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[2, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.10 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(0.10 / model$p[i], times = model$p[i]), phi = rep(0.10 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[3, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.10 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(1.00 / model$p[i], times = model$p[i]), phi = rep(0.10 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[4, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.50 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(0.01 / model$p[i], times = model$p[i]), phi = rep(0.50 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[5, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.50 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(0.10 / model$p[i], times = model$p[i]), phi = rep(0.50 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[6, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.50 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(1.00 / model$p[i], times = model$p[i]), phi = rep(0.50 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[7, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.90 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(0.01 / model$p[i], times = model$p[i]), phi = rep(0.90 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[8, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.90 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(0.10 / model$p[i], times = model$p[i]), phi = rep(0.90 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[9, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i] / (1 - (model$regress == "joint") * 0.90 * (model$q[i] > 0L)), beta = rep(0, times = model$m[i]), alpha = rep(1.00 / model$p[i], times = model$p[i]), phi = rep(0.90 / model$q[i], times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[10, ] <- convert_struc_list_to_theta_vector(lapply(1:info_par$par_num, function (i) { list(omega = par_start[i], beta = rep(0, times = model$m[i]), alpha = rep(0, times = model$p[i]), phi = rep(0, times = model$q[i])) }), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  try_theta[11, ] <- rep(0, times = info_theta$theta_num)
  for (i in 1:try_num) {
    try_theta[i, !is.na(theta_start)] <- theta_start[!is.na(theta_start)]
    try_theta[i, ] <- pmin(try_theta[i, ], pmax(theta_bound_upper - (theta_bound_upper - theta_bound_lower) / 1e3, theta_bound_upper - 1e-9, na.rm = TRUE), na.rm = TRUE)
    try_theta[i, ] <- pmax(try_theta[i, ], pmin(theta_bound_lower + (theta_bound_upper - theta_bound_lower) / 1e3, theta_bound_lower + 1e-9, na.rm = TRUE), na.rm = TRUE)
  }
  try_theta <- try_theta[!duplicated(try_theta), , drop = FALSE]
  try_num <- nrow(try_theta)
  try_obj <- rep(NA_real_, length = try_num)
  est_details <- list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, info_theta = info_theta, print_progress = print_progress)
  for (i in 1:try_num) {
    try_obj[i] <- likelihood_objective(theta = try_theta[i, ], est_details = est_details)
  }
  if (any(try_obj < 1e100)) {
    status_start <- "success"
    theta_start <- try_theta[which.min(try_obj), ]
  } else if (any(!is.na(try_obj))) {
    status_start <- "zero_likelihood"
    theta_start <- try_theta[which.min(try_obj), ]
  } else {
    status_start <- "failure"
    theta_start <- rep(NA_real_, times = info_theta$theta_num)
  }
  return(list(status_start = status_start, theta_start = theta_start))
}
# ------------------------------------------------------------------------------


