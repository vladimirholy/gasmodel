
# OPTIMIZATION AND HESSIAN WRAPPER FUNCTIONS

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
#' @param ... Additional arguments to be passed to the optimization function. These include arguments required by the objective function, namely \code{est_details} and \code{print_progress}.
#'
#' @return A list with components:
#'   \item{status_optim}{The status of the optimization computation.}
#'   \item{theta_optim}{The optimal solution.}
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
#' \code{\link[gasmodel:wrappers_hessian]{wrappers_hessian}}
#'
#' @name wrappers_optim
NULL
# ------------------------------------------------------------------------------


# Wrappers for Optimization Functions ------------------------------------------
#' @describeIn wrappers_optim Wrapper for optimization function \code{\link[stats:optim]{stats::optim()}}.
#' @export
wrapper_optim_stats <- function(obj_fun, theta_start, theta_bound_lower, theta_bound_upper, ...) {
  if (all(theta_bound_lower == -Inf) && all(theta_bound_upper == Inf)) {
    optim_res <- stats::optim(par = theta_start, fn = obj_fun, ...)
  } else {
    optim_res <- stats::optim(par = theta_start, fn = obj_fun, lower = theta_bound_lower, upper = theta_bound_upper, ...)
  }
  if (inherits(optim_res, 'try-error')) {
    status_hessian <- "failure"
    theta_hessian <- matrix(NA_real_, nrow = length(theta_start))
  } else {
    status_optim <- switch(as.character(optim_res$convergence), "0" = "success", "1" = "iteration_limit_reached",  "10" = "degeneracy_reached", "unknown_status")
    theta_optim <- optim_res$par
  }
  return(list(status_optim = status_optim, theta_optim = theta_optim))
}
# ------------------------------------------------------------------------------


# Wrap the Optimization Function from nloptr Package ---------------------------
#' @describeIn wrappers_optim Wrapper for optimization function \code{\link[nloptr:nloptr]{nloptr::nloptr()}}.
#' @export
wrapper_optim_nloptr <- function(obj_fun, theta_start, theta_bound_lower, theta_bound_upper, ...) {
  if (all(theta_bound_lower == -Inf) && all(theta_bound_upper == Inf)) {
    optim_res <- nloptr::nloptr(x0 = theta_start, eval_f = obj_fun, ...)
  } else {
    optim_res <- nloptr::nloptr(x0 = theta_start, eval_f = obj_fun, lb = theta_bound_lower, ub = theta_bound_upper, ...)
  }
  if (inherits(optim_res, 'try-error')) {
    status_hessian <- "failure"
    theta_hessian <- matrix(NA_real_, nrow = length(theta_start))
  } else {
    status_optim <- switch(as.character(optim_res$status), "1" = "success", "2" = "desired_objective_reached",  "3" = "objective_tolerance_reached", "4" = "variables_tolerance_reached", "5" = "iteration_limit_reached", "6" = "maximum_time_reached", "unknown_status")
    theta_optim <- optim_res$solution
  }
  return(list(status_optim = status_optim, theta_optim = theta_optim))
}
# ------------------------------------------------------------------------------


# Documentation for All Hessian Wrappers ---------------------------------------
#' @title
#' Wrappers for Hessian Functions
#'
#' @description
#' Wrappers of common R Hessian functions.
#' Their purpose is to be passed as the \code{hessian_function} argument in the \code{\link[gasmodel:gas]{gas()}} function.
#'
#' @param obj_fun An objective function.
#' @param theta_optim A numeric vector of the optimal values of the variables.
#' @param ... Additional arguments to be passed to the Hessian function. These include arguments required by the objective function, namely \code{est_details} and \code{print_progress}.
#'
#' @return A list with components:
#'   \item{status_hessian}{The status of the Hessian computation.}
#'   \item{theta_hessian}{The Hessian matrix.}
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
#' \code{\link[gasmodel:wrappers_optim]{wrappers_optim}}
#'
#' @name wrappers_hessian
NULL
# ------------------------------------------------------------------------------


# Wrap the Hessian Function from stats Package ---------------------------------
#' @describeIn wrappers_hessian Wrapper for Hessian function \code{\link[stats:optimHess]{stats::optimHess()}}.
#' @export
wrapper_hessian_stats <- function(obj_fun, theta_optim, ...) {
  hessian_res <- suppressWarnings(try(stats::optimHess(par = theta_optim, fn = obj_fun, ...), silent = TRUE))
  if (inherits(hessian_res, 'try-error')) {
    status_hessian <- "failure"
    theta_hessian <- matrix(NA_real_, nrow = length(theta_optim), ncol = length(theta_optim))
  } else {
    status_hessian <- "success"
    theta_hessian <- hessian_res
  }
  return(list(status_hessian = status_hessian, theta_hessian = theta_hessian))
}
# ------------------------------------------------------------------------------


# Wrap the Hessian Function from pracma Package --------------------------------
#' @describeIn wrappers_hessian Wrapper for Hessian function \code{\link[pracma:hessian]{pracma::hessian()}}.
#' @export
wrapper_hessian_pracma <- function(obj_fun, theta_optim, ...) {
  hessian_res <- suppressWarnings(try(pracma::hessian(x0 = theta_optim, f = obj_fun, ...), silent = TRUE))
  if (inherits(hessian_res, 'try-error')) {
    status_hessian <- "failure"
    theta_hessian <- matrix(NA_real_, nrow = length(theta_optim), ncol = length(theta_optim))
  } else {
    status_hessian <- "success"
    theta_hessian <- hessian_res
  }
  return(list(status_hessian = status_hessian, theta_hessian = theta_hessian))
}
# ------------------------------------------------------------------------------


# Wrap the Hessian Function from numDeriv Package ------------------------------
#' @describeIn wrappers_hessian Wrapper for Hessian function \code{\link[numDeriv:hessian]{numDeriv::hessian()}}.
#' @export
wrapper_hessian_numderiv <- function(obj_fun, theta_optim, ...) {
  hessian_res <- suppressWarnings(try(numDeriv::hessian(x = theta_optim, func = obj_fun, ...), silent = TRUE))
  if (inherits(hessian_res, 'try-error')) {
    status_hessian <- "failure"
    theta_hessian <- matrix(NA_real_, nrow = length(theta_optim), ncol = length(theta_optim))
  } else {
    status_hessian <- "success"
    theta_hessian <- hessian_res
  }
  return(list(status_hessian = status_hessian, theta_hessian = theta_hessian))
}
# ------------------------------------------------------------------------------


