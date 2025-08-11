
# FUNCTIONS FOR BOOTSTRAPPING


# Parametric Bootstrap for Static Model ----------------------------------------
boot_param_static <- function(id, run_details) {
  data <- run_details$data
  model <- run_details$model
  control <- run_details$control
  fun <- run_details$fun
  info_distr <- run_details$info_distr
  info_par <- run_details$info_par
  info_coef <- run_details$info_coef
  comp <- run_details$comp
  pre_num <- max(c(model$p, model$q, 1L))
  burn_num <- 10L + max(model$p) + max(model$q)
  full_num <- pre_num + burn_num + model$t
  average_x <- lapply(1:info_par$par_num, function(i) { colMeans(data$x[[i]], na.rm = TRUE) })
  y <- matrix(NA_real_, nrow = full_num, ncol = model$n)
  x <- lapply(1:info_par$par_num, function(i) { rbind(matrix(NA_real_, nrow = pre_num, ncol = model$m[i]), matrix(average_x[[i]], nrow = burn_num, ncol = model$m[i], byrow = TRUE), data$x[[i]]) })
  par_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  score_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  idx_na <- which(rowSums(sapply(x, function(e) { rowSums(is.na(e)) })) > 0L)
  idx_ok <- which(!((1:full_num) %in% idx_na))
  struc <- convert_coef_vector_to_struc_list(coef_vec = model$coef_est, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  omega_vector <- sapply(struc, function(e) { e$omega })
  beta_list <- lapply(struc, function(e) { e$beta })
  alpha_list <- lapply(struc, function(e) { e$alpha })
  phi_list <- lapply(struc, function(e) { e$phi })
  par_init <- model$par_init
  if (any(is.na(par_init))) {
    par_unc <- sapply(1:info_par$par_num, function(i) { (omega_vector[i] + average_x[[i]] %*% beta_list[[i]]) })
    par_init[is.na(par_init)] <- par_unc[is.na(par_init)]
  }
  if (length(idx_na) > 0L) {
    par_tv[idx_na, ] <- matrix(par_init, nrow = length(idx_na), ncol = info_par$par_num, byrow = TRUE)
    score_tv[idx_na, ] <- 0
  }
  par_tv[idx_ok, ] <- matrix(omega_vector, nrow = length(idx_ok), ncol = info_par$par_num, byrow = TRUE)
  if (any(model$m > 0L)) {
    par_tv[idx_ok, ] <- par_tv[idx_ok, ] + sapply(1L:info_par$par_num, function(i) { x[[i]][idx_ok, , drop = FALSE] %*% beta_list[[i]] })
  }
  for (j in idx_ok) {
    y[j, ] <- fun$random(t = 1L, f = par_tv[j, , drop = FALSE])
    score_tv[j, ] <- fun$score(y = y[j, , drop = FALSE], f = par_tv[j, , drop = FALSE])
  }
  comp$est_details$data$y <- y[(pre_num + burn_num + 1L):full_num, , drop = FALSE]
  result_optim <- be_silent_but_theta(do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments)))
  coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Parametric Bootstrap for Dynamic Joint Model ---------------------------------
boot_param_joint <- function(id, run_details) {
  data <- run_details$data
  model <- run_details$model
  control <- run_details$control
  fun <- run_details$fun
  info_distr <- run_details$info_distr
  info_par <- run_details$info_par
  info_coef <- run_details$info_coef
  comp <- run_details$comp
  pre_num <- max(c(model$p, model$q, 1L))
  burn_num <- 10L + max(model$p) + max(model$q)
  full_num <- pre_num + burn_num + model$t
  average_x <- lapply(1:info_par$par_num, function(i) { colMeans(data$x[[i]], na.rm = TRUE) })
  y <- matrix(NA_real_, nrow = full_num, ncol = model$n)
  x <- lapply(1:info_par$par_num, function(i) { rbind(matrix(NA_real_, nrow = pre_num, ncol = model$m[i]), matrix(average_x[[i]], nrow = burn_num, ncol = model$m[i], byrow = TRUE), data$x[[i]]) })
  par_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  score_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  idx_na <- which(rowSums(sapply(x, function(e) { rowSums(is.na(e)) })) > 0L)
  idx_ok <- which(!((1:full_num) %in% idx_na))
  struc <- convert_coef_vector_to_struc_list(coef_vec = model$coef_est, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  omega_vector <- sapply(struc, function(e) { e$omega })
  beta_list <- lapply(struc, function(e) { e$beta })
  alpha_list <- lapply(struc, function(e) { e$alpha })
  phi_list <- lapply(struc, function(e) { e$phi })
  par_init <- model$par_init
  if (any(is.na(par_init))) {
    par_unc <- sapply(1:info_par$par_num, function(i) { (omega_vector[i] + average_x[[i]] %*% beta_list[[i]]) / (1 - sum(phi_list[[i]])) })
    par_init[is.na(par_init)] <- par_unc[is.na(par_init)]
  }
  if (length(idx_na) > 0L) {
    par_tv[idx_na, ] <- matrix(par_init, nrow = length(idx_na), ncol = info_par$par_num, byrow = TRUE)
    score_tv[idx_na, ] <- 0
  }
  par_tv[idx_ok, ] <- matrix(omega_vector, nrow = length(idx_ok), ncol = info_par$par_num, byrow = TRUE)
  if (any(model$m > 0L)) {
    par_tv[idx_ok, ] <- par_tv[idx_ok, ] + sapply(1L:info_par$par_num, function(i) { x[[i]][idx_ok, , drop = FALSE] %*% beta_list[[i]] })
  }
  cur_e <- rep(NA_real_, info_par$par_num)
  for (j in idx_ok) {
    for (i in 1:info_par$par_num) {
      cur_e[i] <- sum(par_tv[j - seq_along(phi_list[[i]]), i] * phi_list[[i]]) + sum(score_tv[j - seq_along(alpha_list[[i]]), i] * alpha_list[[i]])
    }
    par_tv[j, ] <- par_tv[j, ] + cur_e
    y[j, ] <- fun$random(t = 1L, f = par_tv[j, , drop = FALSE])
    score_tv[j, ] <- fun$score(y = y[j, , drop = FALSE], f = par_tv[j, , drop = FALSE])
  }
  comp$est_details$data$y <- y[(pre_num + burn_num + 1L):full_num, , drop = FALSE]
  result_optim <- be_silent_but_theta(do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments)))
  coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Parametric Bootstrap for Dynamic Separate Model ------------------------------
boot_param_sep <- function(id, run_details) {
  data <- run_details$data
  model <- run_details$model
  control <- run_details$control
  fun <- run_details$fun
  info_distr <- run_details$info_distr
  info_par <- run_details$info_par
  info_coef <- run_details$info_coef
  comp <- run_details$comp
  pre_num <- max(c(model$p, model$q, 1L))
  burn_num <- 10L + max(model$p) + max(model$q)
  full_num <- pre_num + burn_num + model$t
  average_x <- lapply(1:info_par$par_num, function(i) { colMeans(data$x[[i]], na.rm = TRUE) })
  y <- matrix(NA_real_, nrow = full_num, ncol = model$n)
  x <- lapply(1:info_par$par_num, function(i) { rbind(matrix(NA_real_, nrow = pre_num, ncol = model$m[i]), matrix(average_x[[i]], nrow = burn_num, ncol = model$m[i], byrow = TRUE), data$x[[i]]) })
  par_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  score_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  idx_na <- which(rowSums(sapply(x, function(e) { rowSums(is.na(e)) })) > 0L)
  idx_ok <- which(!((1:full_num) %in% idx_na))
  struc <- convert_coef_vector_to_struc_list(coef_vec = model$coef_est, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  omega_vector <- sapply(struc, function(e) { e$omega })
  beta_list <- lapply(struc, function(e) { e$beta })
  alpha_list <- lapply(struc, function(e) { e$alpha })
  phi_list <- lapply(struc, function(e) { e$phi })
  err_tv <- matrix(NA_real_, nrow = full_num, ncol = info_par$par_num)
  err_init <- model$par_init
  par_init <- model$par_init
  if (any(is.na(model$par_init))) {
    par_unc <- sapply(1:info_par$par_num, function(i) { (omega_vector[i] + average_x[[i]] %*% beta_list[[i]]) })
    par_init[is.na(par_init)] <- par_unc[is.na(par_init)]
    err_init[is.na(err_init)] <- 0
  }
  if (length(idx_na) > 0L) {
    err_tv[idx_na, ] <- matrix(err_init, nrow = length(idx_na), ncol = info_par$par_num, byrow = TRUE)
    par_tv[idx_na, ] <- matrix(par_init, nrow = length(idx_na), ncol = info_par$par_num, byrow = TRUE)
    score_tv[idx_na, ] <- 0
  }
  par_tv[idx_ok, ] <- matrix(omega_vector, nrow = length(idx_ok), ncol = info_par$par_num, byrow = TRUE)
  if (any(model$m > 0L)) {
    par_tv[idx_ok, ] <- par_tv[idx_ok, ] + sapply(1L:info_par$par_num, function(i) { x[[i]][idx_ok, , drop = FALSE] %*% beta_list[[i]] })
  }
  cur_e <- rep(NA_real_, info_par$par_num)
  for (j in idx_ok) {
    for (i in 1:info_par$par_num) {
      cur_e[i] <- sum(err_tv[j - seq_along(phi_list[[i]]), i] * phi_list[[i]]) + sum(score_tv[j - seq_along(alpha_list[[i]]), i] * alpha_list[[i]])
    }
    err_tv[j, ] <- cur_e
    par_tv[j, ] <- par_tv[j, ] + cur_e
    y[j, ] <- fun$random(t = 1L, f = par_tv[j, , drop = FALSE])
    score_tv[j, ] <- fun$score(y = y[j, , drop = FALSE], f = par_tv[j, , drop = FALSE])
  }
  comp$est_details$data$y <- y[(pre_num + burn_num + 1L):full_num, , drop = FALSE]
  result_optim <- be_silent_but_theta(do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments)))
  coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Simple Block Boostrap --------------------------------------------------------
boot_block_simple <- function(id, run_details) {
  data <- run_details$data
  model <- run_details$model
  control <- run_details$control
  bootstrap <- run_details$bootstrap
  comp <- run_details$comp
  t <- model$t
  block_length <- bootstrap$block_length
  circle <- c(1:t, 1:t)
  block_num <- round(t / block_length)
  block_range <- rep(block_length, block_num)
  if (block_length * block_num > t) {
    block_range[sample(1:block_num, size = block_length * block_num - t, replace = FALSE)] <- block_length - 1
  } else if (block_length * block_num < t) {
    block_range[sample(1:block_num, size = t - block_length * block_num, replace = FALSE)] <- block_length + 1
  }
  block_start <- (sample(1:block_num, size = block_num, replace = TRUE) - 1) * block_length + 1
  block_end <- block_start + block_range - 1
  idx <- circle[unlist(lapply(1:block_num, FUN = function(i) { block_start[i]:block_end[i] }))]
  comp$est_details$data$y <- data$y[idx, , drop = FALSE]
  comp$est_details$data$x <- lapply(data$x, function(xx) { xx[idx, , drop = FALSE] })
  result_optim <- be_silent_but_theta(do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments)))
  coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Moving Block Boostrap --------------------------------------------------------
boot_block_moving <- function(id, run_details) {
  data <- run_details$data
  model <- run_details$model
  control <- run_details$control
  bootstrap <- run_details$bootstrap
  comp <- run_details$comp
  t <- model$t
  block_length <- bootstrap$block_length
  circle <- c(1:t, 1:t)
  block_num <- round(t / block_length)
  block_range <- rep(block_length, block_num)
  if (block_length * block_num > t) {
    block_range[sample(1:block_num, size = block_length * block_num - t, replace = FALSE)] <- block_length - 1
  } else if (block_length * block_num < t) {
    block_range[sample(1:block_num, size = t - block_length * block_num, replace = FALSE)] <- block_length + 1
  }
  block_start <- sample(1:t, size = block_num, replace = TRUE)
  block_end <- block_start + block_range - 1
  idx <- circle[unlist(lapply(1:block_num, FUN = function(i) { block_start[i]:block_end[i] }))]
  comp$est_details$data$y <- data$y[idx, , drop = FALSE]
  comp$est_details$data$x <- lapply(data$x, function(xx) { xx[idx, , drop = FALSE] })
  result_optim <- be_silent_but_theta(do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments)))
  coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Stationary Block Boostrap ----------------------------------------------------
boot_block_stat <- function(id, run_details) {
  data <- run_details$data
  model <- run_details$model
  control <- run_details$control
  bootstrap <- run_details$bootstrap
  comp <- run_details$comp
  t <- model$t
  block_length <- bootstrap$block_length
  circle <- c(1:t, 1:t)
  block_prob <- 1 / block_length
  block_range <- diff(which(c(TRUE, sample(c(TRUE, FALSE), size = t - 1, replace = TRUE, prob = c(block_prob, 1 - block_prob)), TRUE)))
  block_num <- length(block_range)
  block_start <- sample(1:t, size = block_num, replace = TRUE)
  block_end <- block_start + block_range - 1
  idx <- circle[unlist(lapply(1:block_num, FUN = function(i) { block_start[i]:block_end[i] }))]
  comp$est_details$data$y <- data$y[idx, , drop = FALSE]
  comp$est_details$data$x <- lapply(data$x, function(xx) { xx[idx, , drop = FALSE] })
  result_optim <- be_silent_but_theta(do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments)))
  coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  return(coef_set)
}
# ------------------------------------------------------------------------------


