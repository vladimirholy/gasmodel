
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
#' @param gas_object An optional GAS estimate, i.e. a list of S3 class \code{gas} returned by function \code{\link[gasmodel:gas]{gas()}}.
#' @param method A method used for bootstrapping. Supported methods are \code{"parametric"}, \code{"simple_block"}, \code{"moving_block"}, and \code{"stationary_block"}.
#' @param rep_boot A number of bootstrapping repetitions.
#' @param block_length A length of blocks for methods \code{"simple_block"} and \code{"moving_block"}. A mean length of blocks for method \code{"stationary_block"}.
#' @param quant A numeric vector of probabilities determining quantiles.
#' @param y,x,distr,param,scaling,regress,p,q,par_static,par_link,par_init,lik_skip,coef_fix_value,coef_fix_other,coef_fix_special,coef_bound_lower,coef_bound_upper,coef_est When \code{gas_object} is not supplied, the estimated model can be specified using these individual arguments. See the arguments and value of the \code{\link[gasmodel:gas]{gas()}} function for more details.
#' @param optim_function An optimization function. For suitable wrappers of common R optimization functions, see \code{\link[gasmodel:wrappers_optim]{wrappers_optim}}.
#' @param optim_arguments An optional list of arguments to be passed to the optimization function.
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
#' \donttest{# Load the Daily Toilet Paper Sales dataset
#' data("toilet_paper_sales")
#' y <- toilet_paper_sales$quantity
#' x <- as.matrix(toilet_paper_sales[3:9])
#'
#' # Estimate GAS model based on the negative binomial distribution
#' est_negbin <- gas(y = y, x = x, distr = "negbin", regress = "sep")
#' est_negbin
#'
#' # Bootstrap the model (can be time-consuming for a larger number of samples)
#' boot_negbin <- gas_bootstrap(est_negbin, rep_boot = 10)
#' boot_negbin
#'
#' # Plot boxplot of bootstrapped coefficients
#' plot(boot_negbin)}
#'
#' @export
gas_bootstrap <- function(gas_object = NULL, method = "parametric", rep_boot = 1000L, block_length = NULL, quant = c(0.025, 0.975), y = NULL, x = NULL, distr = NULL, param = NULL, scaling = "unit", regress = "joint", p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, lik_skip = 0L, coef_fix_value = NULL, coef_fix_other = NULL, coef_fix_special = NULL, coef_bound_lower = NULL, coef_bound_upper = NULL, coef_est = NULL, optim_function = wrapper_optim_nloptr, optim_arguments = list(opts = list(algorithm = 'NLOPT_LN_NELDERMEAD', xtol_rel = 0, maxeval = 1e2)), parallel_function = NULL, parallel_arguments = list()) {
  if (!is.null(gas_object) && inherits(gas_object, "gas")) {
    gas_bootstrap(gas_object = NULL, method = method, rep_boot = rep_boot, block_length = block_length, quant = quant, y = gas_object$data$y, x = gas_object$data$x, distr = gas_object$model$distr, param = gas_object$model$param, scaling = gas_object$model$scaling, regress = gas_object$model$regress, p = gas_object$model$p, q = gas_object$model$q, par_static = gas_object$model$par_static, par_link = gas_object$model$par_link, par_init = gas_object$model$par_init, lik_skip = gas_object$model$lik_skip, coef_fix_value = gas_object$model$coef_fix_value, coef_fix_other = gas_object$model$coef_fix_other, coef_fix_special = gas_object$model$coef_fix_special, coef_bound_lower = gas_object$model$coef_bound_lower, coef_bound_upper = gas_object$model$coef_bound_upper, coef_est = gas_object$fit$coef_est, optim_function = gas_object$control$optim_function, optim_arguments = gas_object$control$optim_arguments, parallel_function = parallel_function, parallel_arguments = parallel_arguments)
  } else if (!is.null(gas_object)) {
    stop("Unsupported class of gas_object.")
  } else {
    # Load auxiliary variables:
    load <- load_bootstrap(method = method, rep_boot = rep_boot, block_length = block_length, quant = quant, y = y, x = x, distr = distr, param = param, scaling = scaling, regress = regress, p = p, q = q, par_static = par_static, par_link = par_link, par_init = par_init, lik_skip = lik_skip, coef_fix_value = coef_fix_value, coef_fix_other = coef_fix_other, coef_fix_special = coef_fix_special, coef_bound_lower = coef_bound_lower, coef_bound_upper = coef_bound_upper, coef_est = coef_est, optim_function = optim_function, optim_arguments = optim_arguments, parallel_function = parallel_function, parallel_arguments = parallel_arguments)
    data <- load$data
    model <- load$model
    control <- load$control
    fun <- load$fun
    info_distr <- load$info_distr
    info_par <- load$info_par
    info_coef <- load$info_coef
    info_theta <- load$info_theta
    comp <- load$comp
    bootstrap <- load$bootstrap
    # Set parameteric bootstrap for static model:
    if (bootstrap$method == "parametric" && all(model$p + model$q == 0L)) {
      comp$run_details <- list(data = data, model = model, control = control, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp)
      comp$run_fun <- boot_param_static
    # Set parameteric bootstrap for dynamic joint model:
    } else if (bootstrap$method == "parametric" && model$regress == "joint") {
      comp$run_details <- list(data = data, model = model, control = control, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp)
      comp$run_fun <- boot_param_joint
    # Set parameteric bootstrap for dynamic separate model:
    } else if (bootstrap$method == "parametric" && model$regress == "sep") {
      comp$run_details <- list(data = data, model = model, control = control, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, comp = comp)
      comp$run_fun <- boot_param_sep
    # Set simple block bootstrap:
    } else if (bootstrap$method == "simple_block") {
      bootstrap$block_length <- check_my_block_length(block_length = block_length, t = model$t)
      comp$run_details <- list(data = data, model = model, control = control, bootstrap = bootstrap, comp = comp)
      comp$run_fun <- boot_block_simple
    # Set moving block bootstrap:
    } else if (bootstrap$method == "moving_block") {
      bootstrap$block_length <- check_my_block_length(block_length = block_length, t = model$t)
      comp$run_details <- list(data = data, model = model, control = control, bootstrap = bootstrap, comp = comp)
      comp$run_fun <- boot_block_moving
    # Set stationary block bootstrap:
    } else if (bootstrap$method == "stationary_block") {
      bootstrap$block_length <- check_my_block_length(block_length = block_length, t = model$t)
      comp$run_details <- list(data = data, model = model, control = control, bootstrap = bootstrap, comp = comp)
      comp$run_fun <- boot_block_stat
    }
    # Parallelize computation:
    comp$coef_list <- do.call(control$parallel_function, args = c(list(run_num = comp$rep_boot, run_fun = comp$run_fun, run_details = comp$run_details), control$parallel_arguments))
    # Format results:
    bootstrap$coef_set <- name_matrix(matrix(unlist(comp$coef_list), nrow = comp$rep_boot, byrow = TRUE), paste0("s", 1:comp$rep_boot), info_coef$coef_names, drop = c(FALSE, FALSE))
    info_data <- info_data(y = data$y, x = data$x)
    data$y <- name_matrix(data$y, info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
    data$x <- name_list_of_matrices(data$x, info_par$par_names, info_data$index_time_list, info_data$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
    bootstrap$coef_mean <- name_vector(colMeans(bootstrap$coef_set), info_coef$coef_names)
    bootstrap$coef_vcov <- name_matrix(stats::cov(bootstrap$coef_set), info_coef$coef_names, info_coef$coef_names,  drop = c(FALSE, FALSE))
    bootstrap$coef_sd <-  name_vector(sqrt(diag(bootstrap$coef_vcov)), info_coef$coef_names)
    bootstrap$coef_pval <-  name_vector(apply(bootstrap$coef_set, MARGIN = 2, FUN = function(x) { 1 - 2 * abs(mean(x >= 0) - 0.5) }), info_coef$coef_names)
    bootstrap$coef_quant <- name_matrix(t(matrix(apply(bootstrap$coef_set, 2, stats::quantile, probs = comp$quant), ncol = info_coef$coef_num)), info_coef$coef_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
    # Return results:
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
  be_silent(print(gg_fig))
  invisible(gg_fig)
}
# ------------------------------------------------------------------------------


