
# BOOTSTRAPPING OF GAS MODEL FUNCTION


# Bootstrap GAS Model ----------------------------------------------------------
#' @title Bootstrap GAS Model
#'
#' @description
#' A function for bootstrapping coefficients of generalized autoregressive score (GAS) models of Creal et al. (2013) and Harvey (2013).
#' Method \code{"parametric"} repeatedly simulates time series using the parametric model and re-estimates coefficients.
#' Methods \code{"simple_block"}, \code{"moving_block"}, and \code{"stationary_block"} perform the standard variations of the circular block bootstrap.
#' Instead of supplying arguments about the model, the function can be applied to the \code{gas} object obtained by the \code{\link[gasmodel:gas]{gas()}} function.
#' The function enables parallelization.
#'
#' @inheritParams gas
#' @inheritParams gas_simulate
#'
#' @param method A method used for bootstrapping. Supported methods are \code{"parametric"}, \code{"simple_block"}, \code{"moving_block"}, and \code{"stationary_block"}.
#' @param rep_boot A number of bootstrapping repetitions.
#' @param block_length A length of blocks for methods \code{"simple_block"} and \code{"moving_block"}. A mean length of blocks for method \code{"stationary_block"}.
#' @param quant A numeric vector of probabilities determining quantiles.
#' @param parallel_function A parallelization function. For suitable wrappers of common R parallelization functions, see \code{\link[gasmodel:wrappers_parallel]{wrappers_parallel}}. Can be set to \code{NULL} if no parallelization is to be used.
#' @param parallel_arguments An optional list of arguments to be passed to the optimization function.
#'
#' @return A \code{list} of S3 class \code{gas_bootstrap} with components:
#' \item{data$y}{The time series.}
#' \item{data$x}{The exogenous variables.}
#' \item{model$distr}{The conditional distribution.}
#' \item{model$param}{The parametrization of the conditional distribution.}
#' \item{model$scaling}{The scaling function.}
#' \item{model$regress}{The specification of the regression and dynamic equation.}
#' \item{model$t}{The length of the time series.}
#' \item{model$n}{The dimension of the model.}
#' \item{model$m}{The number of exogenous variables.}
#' \item{model$p}{The score order.}
#' \item{model$q}{The autoregressive order.}
#' \item{model$par_static}{The static parameters.}
#' \item{model$par_link}{The parameters with the logarithmic/logistic links.}
#' \item{model$par_init}{The initial values of the time-varying parameters.}
#' \item{model$lik_skip}{The number of skipped observations at the beginning of the time series or after \code{NA} values in the likelihood computation.}
#' \item{model$coef_fix_value}{The values to which coefficients are fixed.}
#' \item{model$coef_fix_other}{The multiples of the estimated coefficients, which are added to the fixed coefficients.}
#' \item{model$coef_fix_special}{The predefined structures of \code{coef_fix_value} and \code{coef_fix_other}.}
#' \item{model$coef_bound_lower}{The lower bounds on coefficients.}
#' \item{model$coef_bound_upper}{The upper bounds on coefficients.}
#' \item{model$coef_est}{The estimated coefficients.}
#' \item{bootstrap$method}{The method used for bootstrapping.}
#' \item{bootstrap$coef_set}{The bootstrapped sets of coefficients.}
#' \item{bootstrap$coef_mean}{The mean of bootstrapped coefficients.}
#' \item{bootstrap$coef_vcov}{The variance-covariance matrix of bootstrapped coefficients.}
#' \item{bootstrap$coef_sd}{The standard deviation of bootstrapped coefficients.}
#' \item{bootstrap$coef_pval}{The p-value of bootstrapped coefficients.}
#' \item{bootstrap$coef_quant}{The quantiles of bootstrapped coefficients.}
#'
#' @note
#' Supported generic functions for S3 class \code{gas_bootstrap} include \code{\link[base:summary]{summary()}}, \code{\link[base:plot]{plot()}}, \code{\link[stats:coef]{coef()}}, and \code{\link[stats:vcov]{vcov()}}.
#'
#' @references
#' Creal, D., Koopman, S. J., and Lucas, A. (2013). Generalized Autoregressive Score Models with Applications. \emph{Journal of Applied Econometrics}, \strong{28}(5), 777–795. \doi{10.1002/jae.1279}.
#'
#' Harvey, A. C. (2013). \emph{Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series}. Cambridge University Press. \doi{10.1017/cbo9781139540933}.
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
#' \code{\link[gasmodel:wrappers_parallel]{wrappers_parallel}}
#'
#' @examples
#' # Load Level of Lake Huron dataset
#' data(LakeHuron)
#' y <- LakeHuron - 570
#' x <- 1:length(y)
#'
#' # Estimate GAS model based on the normal distribution with dynamic mean
#' est_gas <- gas(y = y, x = x, distr = "norm", regress = "sep",
#'   par_static = c(FALSE, TRUE), coef_start = c(9.99, -0.02, 0.46, 0.67, 0.46))
#' est_gas
#'
#' # Bootstrap the model (can be time-consuming for a larger number of samples)
#' \donttest{boot_gas <- gas_bootstrap(est_gas, rep_boot = 10)
#' boot_gas}
#'
#' # Plot boxplot of bootstrapped coefficients
#' \donttest{plot(boot_gas)}
#'
#' @export
gas_bootstrap <- function(gas_object = NULL, method = "parametric", rep_boot = 1000L, block_length = NULL, quant = c(0.025, 0.975), y = NULL, x = NULL, distr = NULL, param = NULL, scaling = "unit", regress = "joint", p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, lik_skip = 0L, coef_fix_value = NULL, coef_fix_other = NULL, coef_fix_special = NULL, coef_bound_lower = NULL, coef_bound_upper = NULL, coef_est = NULL, optim_function = wrapper_optim_nloptr, optim_arguments = list(opts = list(algorithm = 'NLOPT_LN_NELDERMEAD', xtol_rel = 0, maxeval = 1e6)), parallel_function = NULL, parallel_arguments = list()) {
  if (!is.null(gas_object) && inherits(gas_object, "gas")) {
    gas_bootstrap(gas_object = NULL, method = method, rep_boot = rep_boot, block_length = block_length, quant = quant, y = gas_object$data$y, x = gas_object$data$x, distr = gas_object$model$distr, param = gas_object$model$param, scaling = gas_object$model$scaling, regress = gas_object$model$regress, p = gas_object$model$p, q = gas_object$model$q, par_static = gas_object$model$par_static, par_link = gas_object$model$par_link, par_init = gas_object$model$par_init, lik_skip = gas_object$model$lik_skip, coef_fix_value = gas_object$model$coef_fix_value, coef_fix_other = gas_object$model$coef_fix_other, coef_fix_special = gas_object$model$coef_fix_special, coef_bound_lower = gas_object$model$coef_bound_lower, coef_bound_upper = gas_object$model$coef_bound_upper, coef_est = gas_object$fit$coef_est, optim_function = gas_object$control$optim_function, optim_arguments = gas_object$control$optim_arguments, parallel_function = parallel_function, parallel_arguments = parallel_arguments)
  } else if (!is.null(gas_object)) {
    stop("Unsupported class of gas_object.")
  } else {
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
    if (bootstrap$method == "parametric" && all(model$p + model$q == 0L)) {
      comp$run_details <- list(data = data, model = model, control = control, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp)
      comp$run_fun <- function(id, run_details) {
        data <- run_details$data
        model <- run_details$model
        control <- run_details$control
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
        result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
        coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
        return(coef_set)
      }
    } else if (bootstrap$method == "parametric" && model$regress == "joint") {
      comp$run_details <- list(data = data, model = model, control = control, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp)
      comp$run_fun <- function(id, run_details) {
        data <- run_details$data
        model <- run_details$model
        control <- run_details$control
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
        result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
        coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
        return(coef_set)
      }
    } else if (bootstrap$method == "parametric" && model$regress == "sep") {
      comp$run_details <- list(data = data, model = model, control = control, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp)
      comp$run_fun <- function(id, run_details) {
        data <- run_details$data
        model <- run_details$model
        control <- run_details$control
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
        result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
        coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
        return(coef_set)
      }
    } else if (bootstrap$method == "simple_block") {
      bootstrap$block_length <- check_my_block_length(block_length = block_length, t = model$t)
      comp$run_details <- list(data = data, model = model, control = control, bootstrap = bootstrap, comp = comp)
      comp$run_fun <- function(id, run_details) {
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
        result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
        coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
        return(coef_set)
      }
    } else if (bootstrap$method == "moving_block") {
      bootstrap$block_length <- check_my_block_length(block_length = block_length, t = model$t)
      comp$run_details <- list(data = data, model = model, control = control, bootstrap = bootstrap, comp = comp)
      comp$run_fun <- function(id, run_details) {
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
        result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
        coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
        return(coef_set)
      }
    } else if (bootstrap$method == "stationary_block") {
      bootstrap$block_length <- check_my_block_length(block_length = block_length, t = model$t)
      comp$run_details <- list(data = data, model = model, control = control, bootstrap = bootstrap, comp = comp)
      comp$run_fun <- function(id, run_details) {
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
        result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = comp$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
        coef_set <- convert_theta_vector_to_coef_vector(result_optim$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
        return(coef_set)
      }
    }
    comp$coef_list <- do.call(control$parallel_function, args = c(list(run_num = comp$rep_boot, run_fun = comp$run_fun, run_details = comp$run_details), control$parallel_arguments))
    bootstrap$coef_set <- name_matrix(matrix(unlist(comp$coef_list), nrow = comp$rep_boot, byrow = TRUE), paste0("s", 1:comp$rep_boot), info_coef$coef_names, drop = c(FALSE, FALSE))
    info_data <- info_data(y = data$y, x = data$x)
    data$y <- name_matrix(data$y, info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
    data$x <- name_list_of_matrices(data$x, info_par$par_names, info_data$index_time_list, info_data$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
    bootstrap$coef_mean <- name_vector(colMeans(bootstrap$coef_set), info_coef$coef_names)
    bootstrap$coef_vcov <- name_matrix(stats::cov(bootstrap$coef_set), info_coef$coef_names, info_coef$coef_names,  drop = c(FALSE, FALSE))
    bootstrap$coef_sd <-  name_vector(sqrt(diag(bootstrap$coef_vcov)), info_coef$coef_names)
    bootstrap$coef_pval <-  name_vector(apply(bootstrap$coef_set, MARGIN = 2, FUN = function(x) { 1 - 2 * abs(mean(x >= 0) - 0.5) }), info_coef$coef_names)
    bootstrap$coef_quant <- name_matrix(t(matrix(apply(bootstrap$coef_set, 2, stats::quantile, probs = comp$quant), ncol = info_coef$coef_num)), info_coef$coef_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
    report <- list(data = data, model = model, bootstrap = bootstrap)
    class(report) <- "gas_bootstrap"
    return(report)
  }
}
# ------------------------------------------------------------------------------


# Print Bootstrap Coefficients -------------------------------------------------
#' @export
print.gas_bootstrap <- function(x, ...) {
  info_title <- info_title(distr = x$model$distr, param = x$model$param, scaling = x$model$scaling)
  cat("GAS Model:", info_title$title, "\n")
  cat("\n")
  cat("Method:", switch(x$bootstrap$method, "parametric" = "Parametric Bootstrap"), "\n")
  cat("\n")
  cat("Number of Bootstrap Samples:", nrow(x$bootstrap$coef_set), "\n")
  cat("\n")
  coef_table <- cbind("Original" = x$model$coef_est, "Mean" = x$bootstrap$coef_mean, "Std. Error" = x$bootstrap$coef_sd, "P-Value" = x$bootstrap$coef_pval, x$bootstrap$coef_quant)
  cat("Bootstrapped Coefficients:", "\n")
  print(coef_table)
  invisible(x)
}
# ------------------------------------------------------------------------------


# Summarize Bootstrap Coefficients ---------------------------------------------
#' @export
summary.gas_bootstrap <- function(object, ...) {
  print(object)
  cat("\n")
  cat("Bootstrapped Variance-Covariance Matrix:", "\n")
  print(object$bootstrap$coef_vcov)
  invisible(object)
}
# ------------------------------------------------------------------------------


# Obtain Bootstrapped Coefficients -------------------------------------------------
#' @export
coef.gas_bootstrap <- function(object, ...) {
  coef_set <- object$model$coef_set
  return(coef_set)
}
# ------------------------------------------------------------------------------


# Obtain Bootstrapped Variance-Covariance Matrix -------------------------------
#' @export
vcov.gas_bootstrap <- function(object, ...) {
  coef_vcov <- object$bootstrap$coef_vcov
  return(coef_vcov)
}
# ------------------------------------------------------------------------------


# Plot Bootstrapped Coefficients -----------------------------------------------
#' @importFrom dplyr %>%
#' @importFrom ggplot2 .data
#' @export
plot.gas_bootstrap <- function(x, ...) {
  gg_data <- x$bootstrap$coef_set %>%
    dplyr::as_tibble() %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "coef", values_to = "value") %>%
    dplyr::mutate(coef = factor(.data$coef, levels = unique(.data$coef)))
  gg_fig <- ggplot2::ggplot(gg_data, mapping = ggplot2::aes(.data$coef, .data$value)) +
    ggplot2::geom_boxplot(color = "#800000", fill = "#FFAAAA") +
    ggplot2::labs(title = "Bootstrapped Coefficients", x = "Coefficient", y = "Coefficient Value")
  print(gg_fig)
  invisible(gg_fig)
}
# ------------------------------------------------------------------------------


