
# PARALLELIZATION WRAPPER FUNCTIONS

# Documentation for All Parallelization Wrappers -------------------------------
#' @title
#' Wrappers for Parallelization Functions
#'
#' @description
#' Wrappers of common R parallelization functions.
#' Their purpose is to be passed as the \code{parallel_function} argument in the \code{\link[gasmodel:gas_bootstrap]{gas_bootstrap()}} function.
#'
#' @param run_num A number of iterations.
#' @param run_fun A function to be computed.
#' @param run_details A list of variables used for computation.
#' @param ... Additional arguments to be passed to the parallelization function.
#'
#' @return A list containing computed values.
#'
#' @seealso
#' \code{\link[gasmodel:gas_bootstrap]{gas_bootstrap()}}
#' \code{\link[gasmodel:wrappers_optim]{wrappers_optim}}
#' \code{\link[gasmodel:wrappers_hessian]{wrappers_hessian}}
#'
#' @name wrappers_parallel
NULL
# ------------------------------------------------------------------------------


# Wrap the Standard lapply Function --------------------------------------------
#' @describeIn wrappers_parallel Wrapper for function \code{\link[base:lapply]{base::lapply()}}.
#' @export
wrapper_parallel_none <- function(run_num, run_fun, run_details, ...) {
  run_list <- lapply(X = 1:run_num, FUN = run_fun, run_details = run_details, ...)
  return(run_list)
}
# ------------------------------------------------------------------------------


# Wrap the multicore Parallelization Function from parallel Package ------------
#' @describeIn wrappers_parallel Wrapper for parallelization function \code{\link[parallel:mclapply]{parallel::mclapply()}}.
#' @export
wrapper_parallel_multicore <- function(run_num, run_fun, run_details, ...) {
  run_list <- parallel::mclapply(X = 1:run_num, FUN = run_fun, run_details = run_details, ...)
  return(run_list)
}
# ------------------------------------------------------------------------------


# Wrap the snow Parallelization Function from parallel Package -----------------
#' @describeIn wrappers_parallel Wrapper for parallelization function \code{\link[parallel:mclapply]{parallel::parLapply()}}.
#' @export
wrapper_parallel_snow <- function(run_num, run_fun, run_details, ...) {
  cluster <- parallel::makeCluster(...)
  run_list <- parallel::parLapply(cl = cluster, X = 1:run_num, fun = run_fun, run_details = run_details)
  parallel::stopCluster(cluster)
}
# ------------------------------------------------------------------------------


