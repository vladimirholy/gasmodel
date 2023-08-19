
# OPTIMIZATION WRAPPER FUNCTIONS

# Documentation for All Optimization Wrappers ----------------------------------
#' @title
#' Wrappers for Optimization Functions
#'
#' @description
#' Wrappers of common R optimization functions.
#' Their purpose is to be passed as the \code{optim_function} argument in the \code{\link[gasmodel:gas]{gas()}} function.
#'
#' @param obj_fun An objective function.
#' @param theta_start A numeric vector of starting values of the variables.
#' @param theta_bound_lower A numeric vector of lower bounds on the variables.
#' @param theta_bound_upper A numeric vector of upper bounds on the variables.
#' @param est_details A list of variables used for estimation.
#' @param ... Additional arguments to be passed to the optimization function.
#'
#' @return A list with components:
#'   \item{status_optim}{The status of the optimization computation.}
#'   \item{theta_optim}{The optimal solution.}
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
#' \code{\link[gasmodel:wrappers_hessian]{wrappers_hessian}}
#' \code{\link[gasmodel:wrappers_parallel]{wrappers_parallel}}
#'
#' @name wrappers_optim
NULL
# ------------------------------------------------------------------------------


# Wrap the Optimization Function from stats Package ----------------------------
#' @describeIn wrappers_optim Wrapper for optimization function \code{\link[stats:optim]{stats::optim()}}.
#' @export
wrapper_optim_stats <- function(obj_fun, theta_start, theta_bound_lower, theta_bound_upper, est_details, ...) {
  if (all(theta_bound_lower == -Inf) && all(theta_bound_upper == Inf)) {
    optim_res <- stats::optim(par = theta_start, fn = obj_fun, est_details = est_details, ...)
  } else {
    optim_res <- stats::optim(par = theta_start, fn = obj_fun, lower = theta_bound_lower, upper = theta_bound_upper, est_details = est_details, ...)
  }
  if (inherits(optim_res, 'try-error')) {
    status_optim <- "failure"
    obj_optim <- NA_real_
    theta_optim <- rep(NA_real_, times = length(theta_start))
  } else {
    status_optim <- switch(as.character(optim_res$convergence), "0" = "success", "1" = "iteration_limit_reached",  "10" = "degeneracy_reached", "unknown_status")
    obj_optim <- optim_res$value
    theta_optim <- optim_res$par
  }
  return(list(status_optim = status_optim, obj_optim = obj_optim, theta_optim = theta_optim))
}
# ------------------------------------------------------------------------------


# Wrap the Optimization Function from nloptr Package ---------------------------
#' @describeIn wrappers_optim Wrapper for optimization function \code{\link[nloptr:nloptr]{nloptr::nloptr()}}.
#' @export
wrapper_optim_nloptr <- function(obj_fun, theta_start, theta_bound_lower, theta_bound_upper, est_details, ...) {
  if (all(theta_bound_lower == -Inf) && all(theta_bound_upper == Inf)) {
    optim_res <- nloptr::nloptr(x0 = theta_start, eval_f = obj_fun, est_details = est_details, ...)
  } else {
    optim_res <- nloptr::nloptr(x0 = theta_start, eval_f = obj_fun, lb = theta_bound_lower, ub = theta_bound_upper, est_details = est_details, ...)
  }
  if (inherits(optim_res, 'try-error')) {
    status_optim <- "failure"
    obj_optim <- NA_real_
    theta_optim <- rep(NA_real_, times = length(theta_start))
  } else {
    status_optim <- switch(as.character(optim_res$status), "1" = "success", "2" = "desired_objective_reached",  "3" = "objective_tolerance_reached", "4" = "variables_tolerance_reached", "5" = "iteration_limit_reached", "6" = "maximum_time_reached", "unknown_status")
    obj_optim <- optim_res$objective
    theta_optim <- optim_res$solution
  }
  return(list(status_optim = status_optim, obj_optim = obj_optim, theta_optim = theta_optim))
}
# ------------------------------------------------------------------------------


