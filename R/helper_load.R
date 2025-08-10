
# LOAD VARIABLES FUNCTIONS


# Load Variables for Estimate --------------------------------------------------
load_estimate <- function(y, x, distr, param, scaling, regress, p, q, par_static, par_link, par_init, lik_skip, coef_fix_value, coef_fix_other, coef_fix_special, coef_bound_lower, coef_bound_upper, coef_start, optim_function, optim_arguments, hessian_function, hessian_arguments, print_progress) {
  model <- list()
  model$distr <- check_my_distr(distr = distr)
  model$param <- check_my_param(param = param, distr = model$distr)
  model$scaling <- check_my_scaling(scaling = scaling)
  model$regress <- check_my_regress(regress = regress)
  info_distr <- info_distribution(distr = model$distr, param = model$param)
  data <- list()
  data$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  model$t <- check_my_t(y = data$y)
  model$n <- check_my_n(y = data$y)
  info_par <- info_parameters(distr = model$distr, param = model$param, n = model$n)
  data$x <- check_my_x(x = x, t = model$t, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$m <- check_my_m(x = data$x)
  model$p <- check_my_p(p = p, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$q <- check_my_q(q = q, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static <- check_my_par_static(par_static = par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static[model$m == 0L & model$p == 0L & model$q == 0L] <- TRUE
  data$x[model$par_static] <- list(matrix(NA_real_, nrow = model$t, ncol = 0L))
  model$m[model$par_static] <- 0L
  model$p[model$par_static] <- 0L
  model$q[model$par_static] <- 0L
  model$par_link <- check_my_par_link(par_link = par_link, par_static = model$par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = model$par_link)
  model$m <- name_vector(model$m, info_par$par_names)
  model$p <- name_vector(model$p, info_par$par_names)
  model$q <- name_vector(model$q, info_par$par_names)
  model$par_static <- name_vector(model$par_static, info_par$par_names)
  model$par_link <- name_vector(model$par_link, info_par$par_names)
  model$par_init <- name_vector(check_my_par_init(par_init = par_init, par_num = info_par$par_num), info_par$par_names)
  model$lik_skip <- check_my_lik_skip(lik_skip = lik_skip, t = model$t, p = model$p, q = model$q)
  info_coef <- info_coefficients(m = model$m, p = model$p, q = model$q, par_static = model$par_static, par_names = info_par$par_names, par_num = info_par$par_num, group_names = info_par$group_names, group_of_par_names = info_par$group_of_par_names)
  model$coef_fix_value <- check_my_coef_fix_value(coef_fix_value = coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_other <- check_my_coef_fix_other(coef_fix_other = coef_fix_other, coef_fix_value = model$coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_special <- check_my_coef_fix_special(coef_fix_special = coef_fix_special)
  model[c("coef_fix_value", "coef_fix_other")] <- fixed_coefficients(model = model, info_par = info_par, info_coef = info_coef)
  model$coef_bound_lower <- name_vector(check_my_coef_bound_lower(coef_bound_lower = coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  model$coef_bound_upper <- name_vector(check_my_coef_bound_upper(coef_bound_upper = coef_bound_upper, coef_bound_lower = model$coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  info_theta <- info_thetas(coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other, coef_names = info_coef$coef_names)
  fun <- list()
  fun$loglik <- setup_fun_loglik(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$mean <- setup_fun_mean(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$var <- setup_fun_var(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans, par_static = model$par_static)
  comp <- list()
  comp$coef_start <- check_my_coef_start(coef_start = coef_start, coef_bound_lower = model$coef_bound_lower, coef_bound_upper = model$coef_bound_upper, coef_num = info_coef$coef_num)
  comp$theta_start <- convert_coef_vector_to_theta_vector(comp$coef_start, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$theta_bound_lower <- convert_coef_vector_to_theta_vector(model$coef_bound_lower, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$theta_bound_upper <- convert_coef_vector_to_theta_vector(model$coef_bound_upper, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$compute_start <- any(is.na(comp$theta_start))
  control <- list()
  if (is.null(optim_function)) {
    comp$compute_optim <- FALSE
    control['optim_function'] <- list(NULL)
    control['optim_arguments'] <- list(NULL)
  } else {
    comp$compute_optim <- TRUE
    control$optim_function <- check_generic_function(arg = optim_function, arg_name = "optim_function")
    control$optim_arguments <- check_generic_list(arg = optim_arguments, arg_name = "optim_arguments")
  }
  if (is.null(hessian_function)) {
    comp$compute_hessian <- FALSE
    control['hessian_function'] <- list(NULL)
    control['hessian_arguments'] <- list(NULL)
  } else {
    comp$compute_hessian <- TRUE
    control$hessian_function <- check_generic_function(arg = hessian_function, arg_name = "hessian_function")
    control$hessian_arguments <- check_generic_list(arg = hessian_arguments, arg_name = "hessian_arguments")
  }
  comp$print_progress <- check_generic_logical_scalar(arg = print_progress, arg_name = "print_progress")
  comp$est_details <- list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, print_progress = comp$print_progress)
  return(list(data = data, model = model, control = control, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, info_theta = info_theta, comp = comp))
}
# ------------------------------------------------------------------------------


# Load Variables for Bootstrap -------------------------------------------------
load_bootstrap <- function(method, rep_boot, block_length, quant, y, x, distr, param, scaling, regress, p, q, par_static, par_link, par_init, lik_skip, coef_fix_value, coef_fix_other, coef_fix_special, coef_bound_lower, coef_bound_upper, coef_est, optim_function, optim_arguments, parallel_function, parallel_arguments) {
  model <- list()
  model$distr <- check_my_distr(distr = distr)
  model$param <- check_my_param(param = param, distr = model$distr)
  model$scaling <- check_my_scaling(scaling = scaling)
  model$regress <- check_my_regress(regress = regress)
  info_distr <- info_distribution(distr = model$distr, param = model$param)
  data <- list()
  data$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  model$t <- check_my_t(y = data$y)
  model$n <- check_my_n(y = data$y)
  info_par <- info_parameters(distr = model$distr, param = model$param, n = model$n)
  data$x <- check_my_x(x = x, t = model$t, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$m <- check_my_m(x = data$x)
  model$p <- check_my_p(p = p, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$q <- check_my_q(q = q, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static <- check_my_par_static(par_static = par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static[model$m == 0L & model$p == 0L & model$q == 0L] <- TRUE
  data$x[model$par_static] <- list(matrix(NA_real_, nrow = model$t, ncol = 0L))
  model$m[model$par_static] <- 0L
  model$p[model$par_static] <- 0L
  model$q[model$par_static] <- 0L
  model$par_link <- check_my_par_link(par_link = par_link, par_static = model$par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = model$par_link)
  model$m <- name_vector(model$m, info_par$par_names)
  model$p <- name_vector(model$p, info_par$par_names)
  model$q <- name_vector(model$q, info_par$par_names)
  model$par_static <- name_vector(model$par_static, info_par$par_names)
  model$par_link <- name_vector(model$par_link, info_par$par_names)
  model$par_init <- name_vector(check_my_par_init(par_init = par_init, par_num = info_par$par_num), info_par$par_names)
  model$lik_skip <- check_my_lik_skip(lik_skip = lik_skip, t = model$t, p = model$p, q = model$q)
  info_coef <- info_coefficients(m = model$m, p = model$p, q = model$q, par_static = model$par_static, par_names = info_par$par_names, par_num = info_par$par_num, group_names = info_par$group_names, group_of_par_names = info_par$group_of_par_names)
  model$coef_fix_value <- check_my_coef_fix_value(coef_fix_value = coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_other <- check_my_coef_fix_other(coef_fix_other = coef_fix_other, coef_fix_value = model$coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_special <- check_my_coef_fix_special(coef_fix_special = coef_fix_special)
  model[c("coef_fix_value", "coef_fix_other")] <- fixed_coefficients(model = model, info_par = info_par, info_coef = info_coef)
  model$coef_bound_lower <- name_vector(check_my_coef_bound_lower(coef_bound_lower = coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  model$coef_bound_upper <- name_vector(check_my_coef_bound_upper(coef_bound_upper = coef_bound_upper, coef_bound_lower = model$coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  model$coef_est <- name_vector(check_my_coef_est(coef_est = coef_est, coef_num = info_coef$coef_num), info_coef$coef_names)
  info_theta <- info_thetas(coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other, coef_names = info_coef$coef_names)
  fun <- list()
  fun$loglik <- setup_fun_loglik(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans, par_static = model$par_static)
  fun$random <- setup_fun_random(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  control <- list()
  control$optim_function <- check_generic_function(arg = optim_function, arg_name = "optim_function")
  control$optim_arguments <- check_generic_list(arg = optim_arguments, arg_name = "optim_arguments")
  if (is.null(parallel_function)) {
    control$parallel_function <- wrapper_parallel_none
    control$parallel_arguments <- list()
  } else {
    control$parallel_function <- check_generic_function(arg = parallel_function, arg_name = "parallel_function")
    control$parallel_arguments <- check_generic_list(arg = parallel_arguments, arg_name = "parallel_arguments")
  }
  bootstrap <- list()
  bootstrap$method <- check_my_method(method = method, values = c("parametric", "simple_block", "moving_block", "stationary_block"))
  comp <- list()
  comp$rep_boot <- check_generic_positive_integer_scalar(arg = rep_boot, arg_name = "rep_boot")
  bootstrap$coef_set <- name_matrix(matrix(NA_real_, nrow = comp$rep_boot, ncol = info_coef$coef_num), paste0("coef", 1:comp$rep_boot), info_coef$coef_names)
  comp$quant <- check_generic_probability_vector(arg = quant, arg_name = "quant")
  comp$theta_start <- convert_coef_vector_to_theta_vector(model$coef_est, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$theta_bound_lower <- convert_coef_vector_to_theta_vector(model$coef_bound_lower, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$theta_bound_upper <- convert_coef_vector_to_theta_vector(model$coef_bound_upper, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$est_details <- list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, print_progress = FALSE)
  return(list(data = data, model = model, control = control, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, info_theta = info_theta, comp = comp, bootstrap = bootstrap))
}
# ------------------------------------------------------------------------------


# Load Variables for Filter ----------------------------------------------------
load_filter <- function(method, coef_set, rep_gen, t_ahead, x_ahead, rep_ahead, quant, y, x, distr, param, scaling, regress, p, q, par_static, par_link, par_init, coef_fix_value, coef_fix_other, coef_fix_special, coef_bound_lower, coef_bound_upper, coef_est, coef_vcov) {
  model <- list()
  model$distr <- check_my_distr(distr = distr)
  model$param <- check_my_param(param = param, distr = model$distr)
  model$scaling <- check_my_scaling(scaling = scaling)
  model$regress <- check_my_regress(regress = regress)
  info_distr <- info_distribution(distr = model$distr, param = model$param)
  data <- list()
  data$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  model$t <- check_my_t(y = data$y)
  model$t_ahead <- check_my_t(t = t_ahead, positive = FALSE)
  model$n <- check_my_n(y = data$y)
  info_par <- info_parameters(distr = model$distr, param = model$param, n = model$n)
  data$x <- check_my_x(x = x, t = model$t, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$m <- check_my_m(x = data$x)
  data$x_ahead <- check_my_x(x = x_ahead, t = model$t_ahead, m = model$m, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$p <- check_my_p(p = p, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$q <- check_my_q(q = q, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static <- check_my_par_static(par_static = par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static[model$m == 0L & model$p == 0L & model$q == 0L] <- TRUE
  data$x[model$par_static] <- list(matrix(NA_real_, nrow = model$t, ncol = 0L))
  model$m[model$par_static] <- 0L
  model$p[model$par_static] <- 0L
  model$q[model$par_static] <- 0L
  model$par_link <- check_my_par_link(par_link = par_link, par_static = model$par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = model$par_link)
  model$m <- name_vector(model$m, info_par$par_names)
  model$p <- name_vector(model$p, info_par$par_names)
  model$q <- name_vector(model$q, info_par$par_names)
  model$par_static <- name_vector(model$par_static, info_par$par_names)
  model$par_link <- name_vector(model$par_link, info_par$par_names)
  model$par_init <- name_vector(check_my_par_init(par_init = par_init, par_num = info_par$par_num), info_par$par_names)
  info_coef <- info_coefficients(m = model$m, p = model$p, q = model$q, par_static = model$par_static, par_names = info_par$par_names, par_num = info_par$par_num, group_names = info_par$group_names, group_of_par_names = info_par$group_of_par_names)
  model$coef_fix_value <- check_my_coef_fix_value(coef_fix_value = coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_other <- check_my_coef_fix_other(coef_fix_other = coef_fix_other, coef_fix_value = model$coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_special <- check_my_coef_fix_special(coef_fix_special = coef_fix_special)
  model[c("coef_fix_value", "coef_fix_other")] <- fixed_coefficients(model = model, info_par = info_par, info_coef = info_coef)
  model$coef_bound_lower <- name_vector(check_my_coef_bound_lower(coef_bound_lower = coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  model$coef_bound_upper <- name_vector(check_my_coef_bound_upper(coef_bound_upper = coef_bound_upper, coef_bound_lower = model$coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  info_theta <- info_thetas(coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other, coef_names = info_coef$coef_names)
  fun <- list()
  fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans, par_static = model$par_static)
  fun$random <- setup_fun_random(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  filter <- list()
  filter$method <- check_my_method(method = method, values = c("given_coefs", "simulated_coefs"))
  comp <- list()
  comp$rep_ahead <- check_generic_positive_integer_scalar(arg = rep_ahead, arg_name = "rep_ahead")
  comp$quant <- check_generic_probability_vector(arg = quant, arg_name = "quant")
  return(list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, info_theta = info_theta, comp = comp, filter = filter))
}
# ------------------------------------------------------------------------------


# Load Variables for Forecast --------------------------------------------------
load_forecast <- function(method, t_ahead, x_ahead, rep_ahead, quant, y, x, distr, param, scaling, regress, p, q, par_static, par_link, par_init, coef_est) {
  model <- list()
  model$distr <- check_my_distr(distr = distr)
  model$param <- check_my_param(param = param, distr = model$distr)
  model$scaling <- check_my_scaling(scaling = scaling)
  model$regress <- check_my_regress(regress = regress)
  info_distr <- info_distribution(distr = model$distr, param = model$param)
  data <- list()
  data$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  model$t <- check_my_t(y = data$y)
  model$t_ahead <- check_my_t(t = t_ahead)
  model$n <- check_my_n(y = data$y)
  info_par <- info_parameters(distr = model$distr, param = model$param, n = model$n)
  data$x <- check_my_x(x = x, t = model$t, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$m <- check_my_m(x = data$x)
  data$x_ahead <- check_my_x(x = x_ahead, t = model$t_ahead, m = model$m, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$p <- check_my_p(p = p, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$q <- check_my_q(q = q, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static <- check_my_par_static(par_static = par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static[model$m == 0L & model$p == 0L & model$q == 0L] <- TRUE
  data$x[model$par_static] <- list(matrix(NA_real_, nrow = model$t, ncol = 0L))
  model$m[model$par_static] <- 0L
  model$p[model$par_static] <- 0L
  model$q[model$par_static] <- 0L
  model$par_link <- check_my_par_link(par_link = par_link, par_static = model$par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = model$par_link)
  model$m <- name_vector(model$m, info_par$par_names)
  model$p <- name_vector(model$p, info_par$par_names)
  model$q <- name_vector(model$q, info_par$par_names)
  model$par_static <- name_vector(model$par_static, info_par$par_names)
  model$par_link <- name_vector(model$par_link, info_par$par_names)
  model$par_init <- name_vector(check_my_par_init(par_init = par_init, par_num = info_par$par_num), info_par$par_names)
  info_coef <- info_coefficients(m = model$m, p = model$p, q = model$q, par_static = model$par_static, par_names = info_par$par_names, par_num = info_par$par_num, group_names = info_par$group_names, group_of_par_names = info_par$group_of_par_names)
  model$coef_est <- name_vector(check_my_coef_est(coef_est = coef_est, coef_num = info_coef$coef_num), info_coef$coef_names)
  fun <- list()
  fun$mean <- setup_fun_mean(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans, par_static = model$par_static)
  fun$random <- setup_fun_random(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  comp <- list()
  comp$pre_num <- max(c(model$p, model$q, 1L))
  comp$full_num <- comp$pre_num + model$t + model$t_ahead
  comp$average_x <- lapply(1:info_par$par_num, function(i) { colMeans(data$x[[i]], na.rm = TRUE) })
  comp$y <- rbind(matrix(NA_real_, nrow = comp$pre_num, ncol = model$n), data$y, matrix(NA_real_, nrow = model$t_ahead, ncol = model$n))
  comp$x <- lapply(1:info_par$par_num, function(i) { rbind(matrix(NA_real_, nrow = comp$pre_num, ncol = model$m[i]), data$x[[i]], data$x_ahead[[i]]) })
  comp$par_tv <- matrix(NA_real_, nrow = comp$full_num, ncol = info_par$par_num)
  comp$score_tv <- matrix(NA_real_, nrow = comp$full_num, ncol = info_par$par_num)
  comp$idx_na <- which(rowSums(cbind(c(rowSums(is.na(comp$y))[1L:(comp$pre_num + model$t)], rep(0L, model$t_ahead)), sapply(comp$x, function(e) { rowSums(is.na(e)) }))) > 0L)
  comp$idx_ok <- which(!((1:comp$full_num) %in% comp$idx_na))
  comp$idx_ok_regular <- comp$idx_ok[comp$idx_ok <= comp$pre_num + model$t]
  comp$idx_ok_ahead <- comp$idx_ok[comp$idx_ok > comp$pre_num + model$t]
  comp$struc <- convert_coef_vector_to_struc_list(coef_vec = model$coef_est, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  comp$omega_vector <- sapply(comp$struc, function(e) { e$omega })
  comp$beta_list <- lapply(comp$struc, function(e) { e$beta })
  comp$alpha_list <- lapply(comp$struc, function(e) { e$alpha })
  comp$phi_list <- lapply(comp$struc, function(e) { e$phi })
  forecast <- list()
  forecast$method <- check_my_method(method = method, values = c("mean_path", "simulated_paths"))
  return(list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp, forecast = forecast))
}
# ------------------------------------------------------------------------------


# Load Variables for Simulate --------------------------------------------------
load_simulate <- function(t_sim, x_sim, distr, param, scaling, regress, n, p, q, par_static, par_link, par_init, coef_est) {
  model <- list()
  model$distr <- check_my_distr(distr = distr)
  model$param <- check_my_param(param = param, distr = model$distr)
  model$scaling <- check_my_scaling(scaling = scaling)
  model$regress <- check_my_regress(regress = regress)
  info_distr <- info_distribution(distr = model$distr, param = model$param)
  model$t_sim <- check_my_t(t = t_sim)
  model$n <- check_my_n(n = n, dim = info_distr$dim)
  info_par <- info_parameters(distr = model$distr, param = model$param, n = model$n)
  data <- list()
  data$x_sim <- check_my_x(x = x_sim, t = model$t_sim, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$m <- check_my_m(x = data$x_sim)
  model$p <- check_my_p(p = p, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$q <- check_my_q(q = q, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static <- check_my_par_static(par_static = par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$par_static[model$m == 0L & model$p == 0L & model$q == 0L] <- TRUE
  data$x_sim[model$par_static] <- list(matrix(NA_real_, nrow = model$t, ncol = 0L))
  model$m[model$par_static] <- 0L
  model$p[model$par_static] <- 0L
  model$q[model$par_static] <- 0L
  model$par_link <- check_my_par_link(par_link = par_link, par_static = model$par_static, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = model$par_link)
  model$m <- name_vector(model$m, info_par$par_names)
  model$p <- name_vector(model$p, info_par$par_names)
  model$q <- name_vector(model$q, info_par$par_names)
  model$par_static <- name_vector(model$par_static, info_par$par_names)
  model$par_link <- name_vector(model$par_link, info_par$par_names)
  model$par_init <- name_vector(check_my_par_init(par_init = par_init, par_num = info_par$par_num), info_par$par_names)
  info_coef <- info_coefficients(m = model$m, p = model$p, q = model$q, par_static = model$par_static, par_names = info_par$par_names, par_num = info_par$par_num, group_names = info_par$group_names, group_of_par_names = info_par$group_of_par_names)
  model$coef_est <- name_vector(check_my_coef_est(coef_est = coef_est, coef_num = info_coef$coef_num), info_coef$coef_names)
  fun <- list()
  fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans, par_static = model$par_static)
  fun$random <- setup_fun_random(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  comp <- list()
  comp$pre_num <- max(c(model$p, model$q, 1L))
  comp$burn_num <- 10L + max(model$p) + max(model$q)
  comp$full_num <- comp$pre_num + comp$burn_num + model$t_sim
  comp$average_x <- lapply(1:info_par$par_num, function(i) { colMeans(data$x_sim[[i]], na.rm = TRUE) })
  comp$y <- matrix(NA_real_, nrow = comp$full_num, ncol = model$n)
  comp$x <- lapply(1:info_par$par_num, function(i) { rbind(matrix(NA_real_, nrow = comp$pre_num, ncol = model$m[i]), matrix(comp$average_x[[i]], nrow = comp$burn_num, ncol = model$m[i], byrow = TRUE), data$x_sim[[i]]) })
  comp$par_tv <- matrix(NA_real_, nrow = comp$full_num, ncol = info_par$par_num)
  comp$score_tv <- matrix(NA_real_, nrow = comp$full_num, ncol = info_par$par_num)
  comp$idx_na <- which(rowSums(sapply(comp$x, function(e) { rowSums(is.na(e)) })) > 0L | c(rep(TRUE, times = comp$pre_num), rep(FALSE, times = comp$burn_num + model$t_sim)))
  comp$idx_ok <- which(!((1:comp$full_num) %in% comp$idx_na))
  comp$struc <- convert_coef_vector_to_struc_list(coef_vec = model$coef_est, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  comp$omega_vector <- sapply(comp$struc, function(e) { e$omega })
  comp$beta_list <- lapply(comp$struc, function(e) { e$beta })
  comp$alpha_list <- lapply(comp$struc, function(e) { e$alpha })
  comp$phi_list <- lapply(comp$struc, function(e) { e$phi })
  return(list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp))
}
# ------------------------------------------------------------------------------


