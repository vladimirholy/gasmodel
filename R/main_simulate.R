
# SIMULATION OF GAS MODEL FUNCTION


# Simulate GAS Model -----------------------------------------------------------
#' @title Simulate GAS Model
#'
#' @description
#' A function for simulation of generalized autoregressive score (GAS) models of Creal et al. (2013) and Harvey (2013).
#' Instead of supplying arguments about the model, the function can be applied to the \code{gas} object obtained by the \code{\link[gasmodel:gas]{gas()}} function.
#'
#' @inheritParams gas
#'
#' @param gas_object An optional GAS estimate, i.e. a list of S3 class \code{gas} returned by function \code{\link[gasmodel:gas]{gas()}}.
#' @param t_sim A number of observations to simulate.
#' @param x_sim Exogenous variables used for simulations. For a single variable common for all time-varying parameters, a numeric vector. For multiple variables common for all time-varying parameters, a numeric matrix with observations in rows. For individual variables for each time-varying parameter, a list of numeric vectors or matrices in the above form. The number of observation must be equal to \code{t_sim}.
#' @param n A dimension of the model. Required only for multivariate models.
#' @param coef_est A numeric vector of estimated coefficients.
#'
#' @return A \code{list} with components:
#' \item{data$x_sim}{The exogenous variables used in simulation.}
#' \item{model$distr}{The conditional distribution.}
#' \item{model$param}{The parametrization of the conditional distribution.}
#' \item{model$scaling}{The scaling function.}
#' \item{model$spec}{The specification of the dynamic equation.}
#' \item{model$t_sim}{The length of the simulated time series.}
#' \item{model$n}{The dimension of the model.}
#' \item{model$m}{The number of exogenous variables.}
#' \item{model$p}{The score order.}
#' \item{model$q}{The autoregressive order.}
#' \item{model$par_static}{The static parameters.}
#' \item{model$par_link}{The parameters with the logarithmic/logistic links.}
#' \item{model$par_init}{The initial values of the time-varying parameters.}
#' \item{model$coef_est}{The estimated coefficients.}
#' \item{simulation$y_sim}{The simulated time series.}
#' \item{simulation$par_tv_sim}{The simulated time-varying parameters.}
#' \item{simulation$score_tv_sim}{The simulated scores.}
#'
#' @references
#' Creal, D., Koopman, S. J., and Lucas, A. (2013). Generalized Autoregressive Score Models with Applications. \emph{Journal of Applied Econometrics}, \strong{28}(5), 777–795. doi: \href{https://doi.org/10.1002/jae.1279}{10.1002/jae.1279}.
#'
#' Harvey, A. C. (2013). \emph{Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series}. Cambridge University Press. doi: \href{https://doi.org/10.1017/cbo9781139540933}{10.1017/cbo9781139540933}.
#'#'
#' @seealso \code{\link[gasmodel:gas]{gas()}}
#'
#' @examples
#' # Simulate GAS model based on the zero-inflated Poisson distr. with dynamic rate
#' sim_zipois <- gas_simulate(t_sim = 50, distr = "zipois",
#'                            par_static = c(FALSE, TRUE),
#'                            coef_est = c(0.2, 0.1, 0.8, 0.2))
#' sim_zipois
#'
#' # Plot the simulated time series
#' plot(sim_zipois$simulation$y_sim, type = "b")
#'
#' @export
gas_simulate <- function(gas_object = NULL, t_sim = 1L, x_sim = NULL, distr = NULL, param = NULL, scaling = "unit", spec = "joint", n = NULL, p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, coef_est = NULL) {
  if (!is.null(gas_object) && "gas" %in% class(gas_object)) {
    gas_simulate(gas_object = NULL, t_sim = t_sim, x_sim = x_sim, distr = gas_object$model$distr, param = gas_object$model$param, scaling = gas_object$model$scaling, spec = gas_object$model$spec, n = gas_object$model$n, p = gas_object$model$p, q = gas_object$model$q, par_static = gas_object$model$par_static, par_link = gas_object$model$par_link, par_init = gas_object$model$par_init, coef_est = gas_object$fit$coef_est)
  } else if (!is.null(gas_object)) {
    stop("Unsupported class of gas_object.")
  } else {
    model <- list()
    model$distr <- check_my_distr(distr = distr)
    model$param <- check_my_param(param = param, distr = model$distr)
    model$scaling <- check_my_scaling(scaling = scaling)
    model$spec <- check_my_spec(spec = spec)
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
    fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans)
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
    if (all(model$p + model$q == 0L)) {
      comp$par_init <- model$par_init
      if (any(is.na(model$par_init))) {
        comp$par_unc <- sapply(1:info_par$par_num, function(i) { (comp$omega_vector[i] + comp$average_x[[i]] %*% comp$beta_list[[i]]) })
        comp$par_init[is.na(comp$par_init)] <- comp$par_unc[is.na(comp$par_init)]
      }
      if (length(comp$idx_na) > 0L) {
        comp$par_tv[comp$idx_na, ] <- matrix(comp$par_init, nrow = length(comp$idx_na), ncol = info_par$par_num, byrow = TRUE)
        comp$score_tv[comp$idx_na, ] <- 0
      }
      comp$par_tv[comp$idx_ok, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok), ncol = info_par$par_num, byrow = TRUE)
      if (any(model$m > 0L)) {
        comp$par_tv[comp$idx_ok, ] <- comp$par_tv[comp$idx_ok, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok, , drop = FALSE] %*% comp$beta_list[[i]] })
      }
      for (j in comp$idx_ok) {
        comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
      }
      for (j in comp$idx_ok) {
        comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
      }
    } else if (model$spec == "joint") {
      comp$par_init <- model$par_init
      if (any(is.na(comp$par_init))) {
        comp$par_unc <- sapply(1:info_par$par_num, function(i) { (comp$omega_vector[i] + comp$average_x[[i]] %*% comp$beta_list[[i]]) / (1 - sum(comp$phi_list[[i]])) })
        comp$par_init[is.na(comp$par_init)] <- comp$par_unc[is.na(comp$par_init)]
      }
      if (length(comp$idx_na) > 0L) {
        comp$par_tv[comp$idx_na, ] <- matrix(comp$par_init, nrow = length(comp$idx_na), ncol = info_par$par_num, byrow = TRUE)
        comp$score_tv[comp$idx_na, ] <- 0
      }
      comp$par_tv[comp$idx_ok, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok), ncol = info_par$par_num, byrow = TRUE)
      if (any(model$m > 0L)) {
        comp$par_tv[comp$idx_ok, ] <- comp$par_tv[comp$idx_ok, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok, , drop = FALSE] %*% comp$beta_list[[i]] })
      }
      cur_e <- rep(NA_real_, info_par$par_num)
      for (j in comp$idx_ok) {
        for (i in 1:info_par$par_num) {
          cur_e[i] <- sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
        }
        comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
        comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
        comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
      }
    } else if (model$spec == "reg_err") {
      comp$err_tv <- matrix(NA_real_, nrow = comp$full_num, ncol = info_par$par_num)
      comp$err_init <- model$par_init
      comp$par_init <- model$par_init
      if (any(is.na(model$par_init))) {
        comp$par_unc <- sapply(1:info_par$par_num, function(i) { (comp$omega_vector[i] + comp$average_x[[i]] %*% comp$beta_list[[i]]) })
        comp$par_init[is.na(comp$par_init)] <- comp$par_unc[is.na(comp$par_init)]
        comp$err_init[is.na(comp$err_init)] <- 0
      }
      if (length(comp$idx_na) > 0L) {
        comp$err_tv[comp$idx_na, ] <- matrix(comp$err_init, nrow = length(comp$idx_na), ncol = info_par$par_num, byrow = TRUE)
        comp$par_tv[comp$idx_na, ] <- matrix(comp$par_init, nrow = length(comp$idx_na), ncol = info_par$par_num, byrow = TRUE)
        comp$score_tv[comp$idx_na, ] <- 0
      }
      comp$par_tv[comp$idx_ok, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok), ncol = info_par$par_num, byrow = TRUE)
      if (any(model$m > 0L)) {
        comp$par_tv[comp$idx_ok, ] <- comp$par_tv[comp$idx_ok, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok, , drop = FALSE] %*% comp$beta_list[[i]] })
      }
      cur_e <- rep(NA_real_, info_par$par_num)
      for (j in comp$idx_ok) {
        for (i in 1:info_par$par_num) {
          cur_e[i] <- sum(comp$err_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
        }
        comp$err_tv[j, ] <- cur_e
        comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
        comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
        comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
      }
    }
    info_data_sim <- info_data(y = matrix(nrow = model$t_sim, ncol = model$n), x = data$x_sim)
    data$x_sim <- name_list_of_matrices(data$x_sim, info_par$par_names, info_data_sim$index_time_list, info_data_sim$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
    simulation <- list()
    simulation$y_sim <- name_matrix(comp$y[(comp$pre_num + comp$burn_num + 1L):comp$full_num, , drop = FALSE], info_data_sim$index_time, info_data_sim$index_series, drop = c(FALSE, TRUE))
    simulation$par_tv_sim <- name_matrix(comp$par_tv[(comp$pre_num + comp$burn_num + 1L):comp$full_num, , drop = FALSE], info_data_sim$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    simulation$score_tv_sim <- name_matrix(comp$score_tv[(comp$pre_num + comp$burn_num + 1L):comp$full_num, , drop = FALSE], info_data_sim$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    report <- list(data = data, model = model, simulation = simulation)
    class(report) <- "gas_simulate"
    return(report)
  }
}
# ------------------------------------------------------------------------------


# Print Simulations ------------------------------------------------------------
#' @export
print.gas_simulate <- function(x, ...) {
  info_title <- info_title(distr = x$model$distr, param = x$model$param, scaling = x$model$scaling)
  cat("GAS Model:", info_title$title, "\n")
  cat("\n")
  cat("Simulations: \n")
  print(x$simulation$y_sim)
  invisible(x)
}
# ------------------------------------------------------------------------------


