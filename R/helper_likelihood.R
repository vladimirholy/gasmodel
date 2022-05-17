
# LIKELIHOOD EVALUATION AND OBJECTIVE FUNCTIONS


# Evaluate Likelihood Function --------------------------------------------------
likelihood_evaluate <- function(coef, data, model, fun, info_distr, info_par, info_coef) {
  pre_num <- max(c(model$p, model$q, 1L))
  data$y <- rbind(matrix(NA, nrow = pre_num, ncol = model$n), data$y)
  data$x <- lapply(1:info_par$par_num, function(i) { rbind(matrix(NA, nrow = pre_num, ncol = model$m[i]), data$x[[i]]) })
  idx_na <- which(rowSums(cbind(rowSums(is.na(data$y)), sapply(data$x, function(e) { rowSums(is.na(e)) }))) > 0L)
  idx_ok <- which(!((1:(model$t + pre_num)) %in% idx_na))
  if (model$lik_skip > 0L) {
    idx_ignr <- sort(unique(unlist(lapply(1L:model$lik_skip, function(i) { c(0, idx_na) + i }))))
    idx_ignr <- idx_ignr[idx_ignr <= (model$t + pre_num) & !(idx_ignr %in% idx_na)]
  } else {
    idx_ignr <- integer(0)
  }
  idx_lik <- (1:(model$t + pre_num))[!(1:(model$t + pre_num) %in% c(idx_na, idx_ignr))]
  struc <- convert_coef_vector_to_struc_list(coef_vec = coef, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  omega_vector <- sapply(struc, function(e) { e$omega })
  beta_list <- lapply(struc, function(e) { e$beta })
  alpha_list <- lapply(struc, function(e) { e$alpha })
  phi_list <- lapply(struc, function(e) { e$phi })
  par_init <- model$par_init
  if (any(is.na(model$par_init))) {
    par_unc <- sapply(1:info_par$par_num, function(i) { (omega_vector[i] + colMeans(data$x[[i]], na.rm = TRUE) %*% beta_list[[i]]) / (1 - sum(phi_list[[i]])) })
    par_init[is.na(par_init)] <- par_unc[is.na(par_init)]
  }
  tv_l <- rep(NA_real_, model$t + pre_num)
  tv_f <- matrix(NA_real_, nrow = model$t + pre_num, ncol = info_par$par_num)
  tv_s <- matrix(NA_real_, nrow = model$t + pre_num, ncol = info_par$par_num)
  if (length(idx_na) > 0L) {
    tv_f[idx_na, ] <- matrix(par_init, nrow = length(idx_na), ncol = info_par$par_num, byrow = TRUE)
    tv_s[idx_na, ] <- 0
  }
  tv_f[idx_ok, ] <- matrix(omega_vector, nrow = length(idx_ok), ncol = info_par$par_num, byrow = TRUE)
  if (any(model$m > 0L)) {
    tv_f[idx_ok, ] <- tv_f[idx_ok, ] + sapply(1L:info_par$par_num, function(i) { data$x[[i]][idx_ok, , drop = FALSE] %*% beta_list[[i]] })
  }
  if (any(model$p + model$q > 0L)) {
    cur_f <- rep(NA_real_, info_par$par_num)
    for (j in idx_ok) {
      for (i in 1:info_par$par_num) {
        cur_f[i] <- tv_f[j, i] + sum(tv_f[j - seq_along(phi_list[[i]]), i] * phi_list[[i]]) + sum(tv_s[j - seq_along(alpha_list[[i]]), i] * alpha_list[[i]])
      }
      tv_f[j, ] <- cur_f
      tv_s[j, ] <- fun$score(y = data$y[j, , drop = FALSE], f = tv_f[j, , drop = FALSE])
    }
  }
  tv_l[idx_lik] <- -Inf
  try(tv_l[idx_lik] <- fun$loglik(y = data$y[idx_lik, , drop = FALSE], f = tv_f[idx_lik, , drop = FALSE]))
  tv_l[is.na(tv_l)] <- -Inf
  tv_l[idx_na] <- NA_real_
  tv_l <- tv_l[-(1:pre_num)]
  tv_f[idx_na, ] <- NA_real_
  tv_f <- tv_f[-(1:pre_num), , drop = FALSE]
  tv_s[idx_na, ] <- NA_real_
  tv_s <- tv_s[-(1:pre_num), , drop = FALSE]
  eval_tv <- list(lik = tv_l, par = tv_f, score = tv_s)
  return(eval_tv)
}
# ------------------------------------------------------------------------------


# Compute Objective Function ---------------------------------------------------
likelihood_objective <- function(theta, est_details, print_progress) {
  coef <- convert_theta_vector_to_coef_vector(theta_vec = theta, coef_fix_value = est_details$model$coef_fix_value, coef_fix_other = est_details$model$coef_fix_other)
  eval_tv <- suppressWarnings(likelihood_evaluate(coef = coef, data = est_details$data, model = est_details$model, fun = est_details$fun, info_distr = est_details$info_distr, info_par = est_details$info_par, info_coef = est_details$info_coef))
  lik_mean <- mean(eval_tv$lik, na.rm = TRUE)
  if (print_progress) {
    message("Theta: ", paste(formatC(theta, digits = 9, format = "f"), collapse = ", "), "; Log-Likelihood: ", formatC(lik_mean, digits = 9, format = "f"))
  }
  obj <- min(-lik_mean, 1e100, na.rm = TRUE)
  return(obj)
}
# ------------------------------------------------------------------------------


