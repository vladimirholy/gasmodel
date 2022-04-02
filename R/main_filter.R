
# FILTERING OF GAS MODEL FUNCTION


# Filter GAS Model -------------------------------------------------------------
#' @title Filter GAS Model
#'
#' @description
#' A function for obtaining filtered time-varying parameters of generalized autoregressive score (GAS) models.
#' It captures parameter uncertainty and can also be used for forecasting.
#' Method \code{"simulated_coefs"} computes a path of time-varying parameters for each simulated coefficient set under assumption of asymptotic normality with given variance-covariance matrix (see Blasques et al., 2016).
#' Method \code{"given_coefs"} computes a path of time-varying parameters for each supplied coefficient set.
#' Instead of supplying arguments about the model, the function can be applied to the \code{gas} object obtained by the \code{\link[gasmodel:gas]{gas()}} function.
#'
#' @inheritParams gas
#' @inheritParams gas_forecast
#'
#' @param method A method used for parameter uncertainty. Supported methods are \code{"given_coefs"} and \code{"simulated_coefs"}.
#' @param coef_set A numeric matrix of coefficient sets in rows for \code{method = "given_coefs"}. Can be generated for example by \code{\link[gasmodel:gas_bootstrap]{gas_bootstrap()}}.
#' @param rep_gen A number of generated coefficient sets for \code{method = "simulated_coefs"}.
#' @param rep_ahead A number of simulation repetitions for forecasting when \code{t_ahead > 0}.
#' @param quant A numeric vector of probabilities determining quantiles.
#' @param coef_vcov A numeric matrix of estimated covariances between coefficients.
#'
#' @return A \code{list} with components:
#' \item{data$y}{The time series.}
#' \item{data$x}{The exogenous variables.}
#' \item{data$x_ahead}{The out-of-sample exogenous variables. Only when \code{t_ahead > 0}.}
#' \item{model$distr}{The conditional distribution.}
#' \item{model$param}{The parametrization of the conditional distribution.}
#' \item{model$scaling}{The scaling function.}
#' \item{model$t}{The length of the time series.}
#' \item{model$t_ahead}{The length of the out-of-sample time series. Only when \code{t_ahead > 0}.}
#' \item{model$n}{The dimension of the model.}
#' \item{model$m}{The number of exogenous variables.}
#' \item{model$p}{The score order.}
#' \item{model$q}{The autoregressive order.}
#' \item{model$par_static}{The static parameters.}
#' \item{model$par_link}{The parameters with the logarithmic/logistic links.}
#' \item{model$par_init}{The initial values of the time-varying parameters.}
#' \item{model$coef_fix_value}{The values to which coefficients are fixed.}
#' \item{model$coef_fix_other}{The multiples of the estimated coefficients, which are added to the fixed coefficients.}
#' \item{model$coef_fix_special}{The predefined structures of \code{coef_fix_value} and \code{coef_fix_other}.}
#' \item{model$coef_bound_lower}{The lower bounds on coefficients.}
#' \item{model$coef_bound_upper}{The upper bounds on coefficients.}
#' \item{model$coef_set}{The coefficient sets.}
#' \item{filter$method}{The method used for parameter uncertainty.}
#' \item{filter$par_tv_mean}{The mean of the time-varying parameters.}
#' \item{filter$par_tv_sd}{The standard deviation of the time-varying parameters.}
#' \item{filter$par_tv_quant}{The quantiles of the time-varying parameters.}
#' \item{filter$score_tv_mean}{The mean of the scores.}
#' \item{filter$score_tv_sd}{The standard deviation of the scores.}
#' \item{filter$score_tv_quant}{The quantiles of the scores.}
#' \item{filter$y_ahead_mean}{The mean of the forecasted time series. Only when \code{t_ahead > 0}.}
#' \item{filter$y_ahead_sd}{The standard deviation of the forecasted time series. Only when \code{t_ahead > 0}.}
#' \item{filter$y_ahead_quant}{The quantiles of the forecasted time series. Only when \code{t_ahead > 0}.}
#' \item{filter$par_tv_ahead_mean}{The mean of the forecasted time-varying parameters. Only when \code{t_ahead > 0}.}
#' \item{filter$par_tv_ahead_sd}{The standard deviation of the forecasted time-varying parameters. Only when \code{t_ahead > 0}.}
#' \item{filter$par_tv_ahead_quant}{The quantiles of the forecasted time-varying parameters. Only when \code{t_ahead > 0}.}
#' \item{filter$score_tv_ahead_mean}{The mean of the forecasted scores. Only when \code{t_ahead > 0}.}
#' \item{filter$score_tv_ahead_sd}{The standard deviation of the forecasted scores. Only when \code{t_ahead > 0}.}
#' \item{filter$score_tv_ahead_quant}{The quantiles of the forecasted scores. Only when \code{t_ahead > 0}.}
#'
#' @references
#' Blasques, F., Koopman, S. J., Łasak, K., and Lucas, A. (2016). In-Sample Confidence Bands and Out-of-Sample Forecast Bands for Time-Varying Parameters in Observation-Driven Models. \emph{International Journal of Forecasting}, \strong{32}(3), 875–887. doi: \href{https://doi.org/10.1016/j.ijforecast.2015.11.018}{10.1016/j.ijforecast.2015.11.018}.
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
#'
#' @examples
#' # Simulate GAS model based on the normal distr. with dynamic volatility
#' norm_sim <- gas_simulate(distr = "norm", t_sim = 100,
#'                          par_static = c(TRUE, FALSE),
#'                          coef_est = c(0.0, 1.0, 0.4, 0.8))
#' norm_sim
#'
#' # Extract the simulated time series
#' y <- norm_sim$simulation$y_sim
#'
#' # Estimate the model
#' norm_est <- gas(distr = "norm", y = y, par_static = c(TRUE, FALSE))
#' norm_est
#'
#' # Filter the time-varying parameters by the "simulated_coefs" method
#' norm_fil <- gas_filter(norm_est)
#' norm_fil
#'
#' @export
gas_filter <- function(gas_object = NULL, method = "simulated_coefs", coef_set = NULL, rep_gen = 1000L, t_ahead = 0L, x_ahead = NULL, rep_ahead = 1000L, quant = c(0.025, 0.975), y = NULL, x = NULL, distr = NULL, param = NULL, scaling = "unit", p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, coef_fix_value = NULL, coef_fix_other = NULL, coef_fix_special = NULL, coef_bound_lower = NULL, coef_bound_upper = NULL, coef_est = NULL, coef_vcov = NULL) {
  if (!is.null(gas_object) && "gas" %in% class(gas_object)) {
    gas_filter(gas_object = NULL, method = method, coef_set = coef_set, rep_gen = rep_gen, t_ahead = t_ahead, x_ahead = x_ahead, rep_ahead = rep_ahead, quant = quant, y = gas_object$data$y, x = gas_object$data$x, distr = gas_object$model$distr, param = gas_object$model$param, scaling = gas_object$model$scaling, p = gas_object$model$p, q = gas_object$model$q, par_static = gas_object$model$par_static, par_link = gas_object$model$par_link, par_init = gas_object$model$par_init, coef_fix_value = gas_object$model$coef_fix_value, coef_fix_other = gas_object$model$coef_fix_other, coef_fix_special = gas_object$model$coef_fix_special, coef_bound_lower = gas_object$model$coef_bound_lower, coef_bound_upper = gas_object$model$coef_bound_upper, coef_est = gas_object$fit$coef_est, coef_vcov = gas_object$fit$coef_vcov)
  } else if (!is.null(gas_object)) {
    stop("Unsupported class of gas_object.")
  } else {
    model <- list()
    model$distr <- check_my_distr(distr = distr)
    model$param <- check_my_param(param = param, distr = model$distr)
    model$scaling <- check_my_scaling(scaling = scaling)
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
    fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans)
    fun$random <- setup_fun_random(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
    filter <- list()
    filter$method <- check_my_method(method = method, values = c("given_coefs", "simulated_coefs"))
    comp <- list()
    comp$rep_ahead <- check_generic_positive_integer_scalar(arg = rep_ahead, arg_name = "rep_ahead")
    comp$quant <- check_generic_probability_vector(arg = quant, arg_name = "quant")
    if (filter$method == "given_coefs") {
      model$coef_set <- name_matrix(check_my_coef_set(coef_set = coef_set, coef_num = info_coef$coef_num), paste0("coef", 1:nrow(coef_set)), info_coef$coef_names)
      comp$rep_gen <- nrow(model$coef_set)
    } else if (method == "simulated_coefs") {
      comp$rep_gen <- check_generic_positive_integer_scalar(arg = rep_gen, arg_name = "rep_gen")
      comp$coef_est <- check_my_coef_est(coef_est = coef_est, coef_num = info_coef$coef_num)
      comp$theta_est <- convert_coef_vector_to_theta_vector(coef_vec = comp$coef_est, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
      comp$coef_vcov <- check_my_coef_vcov(coef_vcov = coef_vcov, coef_num = info_coef$coef_num)
      comp$theta_vcov <- convert_coef_matrix_to_theta_matrix(coef_mat = comp$coef_vcov, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
      comp$theta_bound_lower <- convert_coef_vector_to_theta_vector(model$coef_bound_lower, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
      comp$theta_bound_upper <- convert_coef_vector_to_theta_vector(model$coef_bound_upper, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
      comp$theta_set <- suppressWarnings(mvnfast::rmvn(comp$rep_gen, mu = comp$theta_est, sigma = comp$theta_vcov))
      model$coef_set <- name_matrix(matrix(NA_real_, nrow = comp$rep_gen, ncol = info_coef$coef_num), paste0("coef", 1:comp$rep_gen), info_coef$coef_names)
      for (i in 1:comp$rep_gen) {
        model$coef_set[i, ] <- convert_theta_vector_to_coef_vector(pmax(pmin(comp$theta_set[i, ], comp$theta_bound_upper), comp$theta_bound_lower), coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
      }
    }
    comp$pre_num <- max(c(model$p, model$q))
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
    comp$y_ahead_path <- array(NA_real_, dim = c(comp$rep_gen * comp$rep_ahead, model$t_ahead, model$n))
    comp$par_tv_path <- array(NA_real_, dim = c(comp$rep_gen, model$t, info_par$par_num))
    comp$par_tv_ahead_path <- array(NA_real_, dim = c(comp$rep_gen * comp$rep_ahead, model$t_ahead, info_par$par_num))
    comp$score_tv_path <- array(NA_real_, dim = c(comp$rep_gen, model$t, info_par$par_num))
    comp$score_tv_ahead_path <- array(NA_real_, dim = c(comp$rep_gen * comp$rep_ahead, model$t_ahead, info_par$par_num))
    for (b in 1:comp$rep_gen) {
      comp$struc <- convert_coef_vector_to_struc_list(coef_vec = model$coef_set[b, ], m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
      comp$omega_vector <- sapply(comp$struc, function(e) { e$omega })
      comp$beta_list <- lapply(comp$struc, function(e) { e$beta })
      comp$alpha_list <- lapply(comp$struc, function(e) { e$alpha })
      comp$phi_list <- lapply(comp$struc, function(e) { e$phi })
      comp$par_init <- model$par_init
      if (any(is.na(comp$par_init))) {
        comp$par_unc <- sapply(1:info_par$par_num, function(i) { (comp$omega_vector[i] + comp$average_x[[i]] %*% comp$beta_list[[i]]) / (1 - sum(comp$phi_list[[i]])) })
        comp$par_init[is.na(comp$par_init)] <- comp$par_unc[is.na(comp$par_init)]
      }
      if (length(comp$idx_na) > 0L) {
        comp$par_tv[comp$idx_na, ] <- matrix(comp$par_init, nrow = length(comp$idx_na), ncol = info_par$par_num, byrow = TRUE)
        comp$score_tv[comp$idx_na, ] <- 0
      }
      comp$par_tv[comp$idx_ok_regular, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok_regular), ncol = info_par$par_num, byrow = TRUE)
      if (any(model$m > 0L)) {
        comp$par_tv[comp$idx_ok_regular, ] <- comp$par_tv[comp$idx_ok_regular, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_regular, , drop = FALSE] %*% comp$beta_list[[i]] })
      }
      cur_f <- rep(NA_real_, info_par$par_num)
      for (j in comp$idx_ok_regular) {
        for (i in 1:info_par$par_num) {
          cur_f[i] <- comp$par_tv[j, i] + sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
        }
        comp$par_tv[j, ] <- cur_f
        comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
      }
      if (model$t_ahead > 0L) {
        for (a in 1:comp$rep_ahead) {
          comp$par_tv[comp$idx_ok_ahead, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok_ahead), ncol = info_par$par_num, byrow = TRUE)
          if (any(model$m > 0L)) {
            comp$par_tv[comp$idx_ok_ahead, ] <- comp$par_tv[comp$idx_ok_ahead, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_ahead, , drop = FALSE] %*% comp$beta_list[[i]] })
          }
          cur_f <- rep(NA_real_, info_par$par_num)
          for (j in comp$idx_ok_ahead) {
            for (i in 1:info_par$par_num) {
              cur_f[i] <- comp$par_tv[j, i] + sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
            }
            comp$par_tv[j, ] <- cur_f
            comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
            comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
          }
          comp$y_ahead_path[(b - 1L) * comp$rep_ahead + a, , ] <- comp$y[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$par_tv_ahead_path[(b - 1L) * comp$rep_ahead + a, , ] <- comp$par_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$score_tv_ahead_path[(b - 1L) * comp$rep_ahead + a, , ] <- comp$score_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
        }
      }
      comp$par_tv_path[b, , ] <- comp$par_tv[(comp$pre_num + 1L):(comp$pre_num + model$t), ]
      comp$score_tv_path[b, , ] <- comp$score_tv[(comp$pre_num + 1L):(comp$pre_num + model$t), ]
    }
    info_data <- info_data(y = data$y, x = data$x)
    data$y <- name_matrix(data$y, info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
    data$x <- name_list_of_matrices(data$x, info_par$par_names, info_data$index_time_list, info_data$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
    filter$par_tv_mean <- name_matrix(colMeans(comp$par_tv_path, na.rm = TRUE), info_data$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    filter$par_tv_sd <- name_matrix(apply(comp$par_tv_path, 2:3, stats::sd, na.rm = TRUE), info_data$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    filter$par_tv_quant <- name_3d_array(aperm(array(apply(comp$par_tv_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t, info_par$par_num)), c(2, 3, 1)), info_data$index_time, info_par$par_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
    filter$score_tv_mean <- name_matrix(colMeans(comp$score_tv_path, na.rm = TRUE), info_data$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    filter$score_tv_sd <- name_matrix(apply(comp$score_tv_path, 2:3, stats::sd, na.rm = TRUE), info_data$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    filter$score_tv_quant <- name_3d_array(aperm(array(apply(comp$score_tv_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t, info_par$par_num)), c(2, 3, 1)), info_data$index_time, info_par$par_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
    if (model$t_ahead > 0L) {
      info_data_ahead <- info_data(y = matrix(nrow = model$t_ahead, ncol = model$n), x = data$x_ahead, skip_t = model$t)
      data$x_ahead <- name_list_of_matrices(data$x_ahead, info_par$par_names, info_data_ahead$index_time_list, info_data_ahead$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
      filter$y_ahead_mean <- name_matrix(colMeans(comp$y_ahead_path, na.rm = TRUE), info_data_ahead$index_time, info_data_ahead$index_series, drop = c(FALSE, TRUE))
      filter$y_ahead_sd <- name_matrix(apply(comp$y_ahead_path, 2:3, stats::sd, na.rm = TRUE), info_data_ahead$index_time, info_data_ahead$index_series, drop = c(FALSE, TRUE))
      filter$y_ahead_quant <- name_3d_array(aperm(array(apply(comp$y_ahead_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t_ahead, model$n)), c(2, 3, 1)), info_data_ahead$index_time, info_data_ahead$index_series, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
      filter$par_tv_ahead_mean <- name_matrix(colMeans(comp$par_tv_ahead_path, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      filter$par_tv_ahead_sd <- name_matrix(apply(comp$par_tv_ahead_path, 2:3, stats::sd, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      filter$par_tv_ahead_quant <- name_3d_array(aperm(array(apply(comp$par_tv_ahead_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t_ahead, info_par$par_num)), c(2, 3, 1)), info_data_ahead$index_time, info_par$par_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
      filter$score_tv_ahead_mean <- name_matrix(colMeans(comp$score_tv_ahead_path, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      filter$score_tv_ahead_sd <- name_matrix(apply(comp$score_tv_ahead_path, 2:3, stats::sd, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      filter$score_tv_ahead_quant <- name_3d_array(aperm(array(apply(comp$score_tv_ahead_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t_ahead, info_par$par_num)), c(2, 3, 1)), info_data_ahead$index_time, info_par$par_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
    } else {
      data$x_ahead <- NULL
      model$t_ahead <- NULL
    }
    report <- list(data = data, model = model, filter = filter)
    class(report) <- "gas_filter"
    return(report)
  }
}
# ------------------------------------------------------------------------------


# Print Filter -----------------------------------------------------------------
#' @export
print.gas_filter <- function(x, ...) {
  info_title <- info_title(distr = x$model$distr, param = x$model$param, scaling = x$model$scaling)
  cat("GAS Model:", info_title$title, "\n")
  cat("\n")
  cat("Method:", switch(x$filter$method, "given_coefs" = "Simulated Paths from Given Coefficients", "simulated_coefs" = "Simulated Paths from Simulated Coefficients"), "\n")
  cat("\n")
  cat("Filtered Parameters: \n")
  print(abind::abind(x$filter$par_tv_quant, x$filter$par_tv_ahead_quant, along = 1L))
  invisible(x)
}
# ------------------------------------------------------------------------------


