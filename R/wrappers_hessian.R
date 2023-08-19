
# HESSIAN WRAPPER FUNCTIONS

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
#' @param est_details A list of variables used for estimation.
#' @param ... Additional arguments to be passed to the Hessian function.
#'
#' @return A list with components:
#'   \item{status_hessian}{The status of the Hessian computation.}
#'   \item{theta_hessian}{The Hessian matrix.}
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
#' \code{\link[gasmodel:wrappers_optim]{wrappers_optim}}
#' \code{\link[gasmodel:wrappers_parallel]{wrappers_parallel}}
#'
#' @name wrappers_hessian
NULL
# ------------------------------------------------------------------------------


# Wrap the Hessian Function from stats Package ---------------------------------
#' @describeIn wrappers_hessian Wrapper for Hessian function \code{\link[stats:optimHess]{stats::optimHess()}}.
#' @export
wrapper_hessian_stats <- function(obj_fun, theta_optim, est_details, ...) {
  hessian_res <- try(stats::optimHess(par = theta_optim, fn = obj_fun, est_details = est_details, ...))
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
wrapper_hessian_pracma <- function(obj_fun, theta_optim, est_details, ...) {
  hessian_res <- try(pracma::hessian(x0 = theta_optim, f = obj_fun, est_details = est_details, ...))
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
wrapper_hessian_numderiv <- function(obj_fun, theta_optim, est_details, ...) {
  hessian_res <- try(numDeriv::hessian(x = theta_optim, func = obj_fun, est_details = est_details, ...))
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


