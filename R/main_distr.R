
# DISTRIBUTION FUNCTIONS


# Get Table of Supported Distributions -----------------------------------------
#' @title Get the Table of Supported Distributions
#'
#' @description
#' A function listing distributions and their parametrizations supported by the \code{\link[gasmodel:gas]{gas()}} function.
#' Output can be filtered using several arguments.
#'
#' @param filter_distr An optional vector of distributions by which the output is filtered.
#' @param filter_param An optional vector of parametrizations by which the output is filtered.
#' @param filter_type An optional vector of data types by which the output is filtered.
#' @param filter_dim An optional vector of dimensions by which the output is filtered.
#' @param filter_orthog An optional logical value indicating whether the parametrization is orthogonal by which the output is filtered.
#' @param filter_default An optional logical value indicating whether the parameterization is the default for the distribution by which the output is filtered.
#'
#' @return A \code{data.frame} with columns:
#' \item{distr_title}{The title of the distribution.}
#' \item{param_title}{The title of the parametrization.}
#' \item{distr}{The distribution.}
#' \item{param}{The parametrization.}
#' \item{type}{The data type.}
#' \item{dim}{The dimension.}
#' \item{orthog}{The indication of whether the parametrization is orthogonal.}
#' \item{default}{The indication of whether the parameterization is the default for the distribution.}
#'
#' @seealso
#' \code{\link[gasmodel:distr_density]{distr_density()}},
#' \code{\link[gasmodel:distr_mean]{distr_mean()}},
#' \code{\link[gasmodel:distr_var]{distr_var()}},
#' \code{\link[gasmodel:distr_score]{distr_score()}},
#' \code{\link[gasmodel:distr_fisher]{distr_fisher()}},
#' \code{\link[gasmodel:distr_random]{distr_random()}},
#' \code{\link[gasmodel:gas]{gas()}}
#'
#' @examples
#' # List all available distributions
#' distr()
#'
#' # List only distributions for count data
#' distr(filter_type = "count")
#'
#' # Show default parametrization for the exponential distribution
#' distr(filter_dist = "exp", filter_default = TRUE)
#'
#' @export
distr <- function(filter_distr = NULL, filter_param = NULL, filter_type = NULL, filter_dim = NULL, filter_orthog = NULL, filter_default = NULL) {
  report_table <- distr_table
  if (is.vector(filter_distr) && !is.list(filter_distr)) {
    report_table <- report_table[report_table$distr %in% filter_distr, ]
  } else if (!is.null(filter_distr)) {
    stop("Ivalid type of argument filter_distr.")
  }
  if (is.vector(filter_param) && !is.list(filter_param)) {
    report_table <- report_table[report_table$param %in% filter_param, ]
  } else if (!is.null(filter_param)) {
    stop("Ivalid type of argument filter_param")
  }
  if (is.vector(filter_type) && !is.list(filter_type)) {
    report_table <- report_table[report_table$type %in% filter_type, ]
  } else if (!is.null(filter_type)) {
    stop("Ivalid type of argument filter_type")
  }
  if (is.vector(filter_dim) && !is.list(filter_dim)) {
    report_table <- report_table[report_table$dim %in% filter_dim, ]
  } else if (!is.null(filter_dim)) {
    stop("Ivalid type of argument filter_dim")
  }
  if (is.vector(filter_orthog) && !is.list(filter_orthog)) {
    report_table <- report_table[report_table$orthog %in% filter_orthog, ]
  } else if (!is.null(filter_orthog)) {
    stop("Ivalid type of argument filter_orthog")
  }
  if (is.vector(filter_default) && !is.list(filter_default)) {
    report_table <- report_table[report_table$default %in% filter_default, ]
  } else if (!is.null(filter_default)) {
    stop("Ivalid type of argument filter_default")
  }
  return(report_table)
}
# ------------------------------------------------------------------------------


# Compute Density --------------------------------------------------------------
#' @title Compute Density
#'
#' @description
#' A function computing density or its logarithm of a given distribution.
#'
#' @param y Observations. For an univariate distribution, a numeric vector. For a multivariate distribution, a numeric matrix with observations in rows or a numeric vector of a single observation.
#' @param f Parameters. For the same parameters for all observations, a numeric vector. For individual parameters for each observation, a numeric matrix with rows corresponding to observations.
#' @param distr A distribution.
#' @param param A parametrization of the distribution.
#' @param par_link An optional logical vector indicating whether the logarithmic/logistic link should be applied to restricted parameters in order to obtain unrestricted values. Defaults to keeping the original link for all parameters.
#' @param trans An optional transformation of the density. The supported transformation is the logarithm of the density (\code{trans = "log"}).
#'
#' @return The (transformed) density.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}}
#'
#' @examples
#' # Density of the negative binomial distribution
#' distr_density(y = c(1, 8, 5, 0, 0), f = c(0.5, 1.2), distr = "negbin")
#'
#' # Density of the multivariate normal distribution
#' distr_density(y = rbind(c(0.5, 0.6), c(-2.3, -1.8), c(-0.2, 0.2)),
#'               f = c(0, 0, 1, 1, 0.5), distr = "mnorm")
#'
#' @export
distr_density <- function(y, f, distr, param = NULL, par_link = NULL, trans = NULL) {
  comp <- list()
  comp$distr <- check_my_distr(distr = distr)
  comp$param <- check_my_param(param = param, distr = comp$distr)
  info_distr <- info_distribution(distr = comp$distr, param = comp$param)
  comp$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  comp$t <- check_my_t(y = comp$y)
  comp$n <- check_my_n(y = comp$y)
  info_par <- info_parameters(distr = comp$distr, param = comp$param, n = comp$n)
  comp$f <- check_my_f(f = f, t = comp$t, par_num = info_par$par_num)
  comp$par_link <- check_my_par_link(par_link = par_link, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = comp$par_link)
  if (is.null(trans)) {
    comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f, par_trans = info_par$par_trans)
    res_density <- do.call(paste("distr", comp$distr, comp$param, "density", sep = "_"), args = list(y = comp$y, f = comp$f_orig))
  } else if (is.vector(trans) && !is.list(trans) && length(trans) == 1L && trans == "log") {
    comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f, par_trans = info_par$par_trans)
    res_density <- do.call(paste("distr", comp$distr, comp$param, "loglik", sep = "_"), args = list(y = comp$y, f = comp$f_orig))
  } else {
    stop("Invalid value of trans.")
  }
  res_density <- as.vector(res_density)
  return(res_density)
}
# ------------------------------------------------------------------------------


# Compute Mean -----------------------------------------------------------------
#' @title Compute Mean
#'
#' @description
#' A function computing mean for a given distribution.
#'
#' @inheritParams distr_density
#'
#' @return The mean.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}}
#'
#' @examples
#' # Mean for the negative binomial distribution
#' distr_mean(f = c(0.5, 1.2), distr = "negbin")
#'
#' # Mean for the multivariate normal distribution
#' distr_mean(f = c(0, 0, 1, 1, 0.5), distr = "mnorm")
#'
#' @export
distr_mean <- function(f, distr, param = NULL, par_link = NULL) {
  comp <- list()
  comp$distr <- check_my_distr(distr = distr)
  comp$param <- check_my_param(param = param, distr = comp$distr)
  info_distr <- info_distribution(distr = comp$distr, param = comp$param)
  comp$f <- check_my_f(f = f)
  comp$t <- check_my_t(y = comp$f)
  comp$n <- check_my_n(f = comp$f, distr = comp$distr, param = comp$param, dim = info_distr$dim)
  info_par <- info_parameters(distr = comp$distr, param = comp$param, n = comp$n)
  comp$par_link <- check_my_par_link(par_link = par_link, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = comp$par_link)
  comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f, par_trans = info_par$par_trans)
  res_mean <- drop(do.call(paste("distr", comp$distr, comp$param, "mean", sep = "_"), args = list(f = comp$f_orig)))
  return(res_mean)
}
# ------------------------------------------------------------------------------


# Compute Variance -------------------------------------------------------------
#' @title Compute Variance
#'
#' @description
#' A function computing variance for a given distribution.
#'
#' @inheritParams distr_density
#'
#' @return The variance.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}}
#'
#' @examples
#' # Variance for the negative binomial distribution
#' distr_var(f = c(0.5, 1.2), distr = "negbin")
#'
#' # Variance for the multivariate normal distribution
#' distr_var(f = c(0, 0, 1, 1, 0.5), distr = "mnorm")
#'
#' @export
distr_var <- function(f, distr, param = NULL, par_link = NULL) {
  comp <- list()
  comp$distr <- check_my_distr(distr = distr)
  comp$param <- check_my_param(param = param, distr = comp$distr)
  info_distr <- info_distribution(distr = comp$distr, param = comp$param)
  comp$f <- check_my_f(f = f)
  comp$t <- check_my_t(y = comp$f)
  comp$n <- check_my_n(f = comp$f, distr = comp$distr, param = comp$param, dim = info_distr$dim)
  info_par <- info_parameters(distr = comp$distr, param = comp$param, n = comp$n)
  comp$par_link <- check_my_par_link(par_link = par_link, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = comp$par_link)
  comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f, par_trans = info_par$par_trans)
  res_var <- drop(do.call(paste("distr", comp$distr, comp$param, "var", sep = "_"), args = list(f = comp$f_orig)))
  return(res_var)
}
# ------------------------------------------------------------------------------


# Compute Score ----------------------------------------------------------------
#' @title Compute Score
#'
#' @description
#' A function computing score or scaled score for a given distribution.
#'
#' @inheritParams distr_density
#'
#' @param scaling An optional scaling function for the score. The supported scaling functions are the unit scaling (\code{scaling = "unit"}), the inverse of the Fisher information matrix scaling (\code{scaling = "fisher_inv"}), and the inverse square root of the Fisher information matrix scaling (\code{scaling = "fisher_inv_sqrt"}).
#'
#' @return The (scaled) score.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}}
#'
#' @examples
#' # Score for the negative binomial distribution
#' distr_score(y = c(1, 8, 5, 0, 0), f = c(0.5, 1.2), distr = "negbin")
#'
#' # Score for the multivariate normal distribution
#' distr_score(y = rbind(c(0.5, 0.6), c(-2.3, -1.8), c(-0.2, 0.2)),
#'             f = c(0, 0, 1, 1, 0.5), distr = "mnorm")
#'
#' @export
distr_score <- function(y, f, distr, param = NULL, par_link = NULL, scaling = NULL) {
  comp <- list()
  comp$distr <- check_my_distr(distr = distr)
  comp$param <- check_my_param(param = param, distr = comp$distr)
  info_distr <- info_distribution(distr = comp$distr, param = comp$param)
  comp$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  comp$t <- check_my_t(y = comp$y)
  comp$n <- check_my_n(y = comp$y)
  info_par <- info_parameters(distr = comp$distr, param = comp$param, n = comp$n)
  comp$f <- check_my_f(f = f, t = comp$t, par_num = info_par$par_num)
  comp$par_link <- check_my_par_link(par_link = par_link, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = comp$par_link)
  comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f, par_trans = info_par$par_trans)
  res_score <- matrix(NA_real_, nrow = comp$t, ncol = info_par$par_num)
  if (is.null(scaling) || (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && scaling == "unit")) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      res_score[i, ] <- as.vector(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "score", sep = "_"), args = list(y = comp$y[i, , drop = FALSE], f = comp$f_orig))[1, ])
    }
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && scaling == "fisher_inv" && info_distr$orthog == TRUE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob <- reparam_link_jacob(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      comp$fisher_inv_orig <- matrix_diag_inv(do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ])
      res_score[i, ] <- as.vector(comp$jacob %*% comp$fisher_inv_orig %*% do.call(paste("distr", comp$distr, comp$param, "score", sep = "_"), args = list(y = comp$y[i, , drop = FALSE], f = comp$f_orig))[1, ])
    }
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && scaling == "fisher_inv" && info_distr$orthog == FALSE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob <- reparam_link_jacob(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      comp$fisher_inv_orig <- matrix_inv(do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ])
      res_score[i, ] <- as.vector(comp$jacob %*% comp$fisher_inv_orig %*% do.call(paste("distr", comp$distr, comp$param, "score", sep = "_"), args = list(y = comp$y[i, , drop = FALSE], f = comp$f_orig))[1, ])
    }
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && scaling == "fisher_inv_sqrt" && info_distr$orthog == TRUE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      comp$fisher_inv_sqrt_tilde <- matrix_diag_inv_sqrt(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv)
      res_score[i, ] <- as.vector(comp$fisher_inv_sqrt_tilde %*% t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "score", sep = "_"), args = list(y = comp$y[i, , drop = FALSE], f = comp$f_orig))[1, ])
    }
  } else if (is.vector(scaling) && !is.list(scaling) && length(scaling) == 1L && scaling == "fisher_inv_sqrt" && info_distr$orthog == FALSE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      comp$fisher_inv_sqrt_tilde <- matrix_inv_sqrt(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv)
      res_score[i, ] <- as.vector(comp$fisher_inv_sqrt_tilde %*% t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "score", sep = "_"), args = list(y = comp$y[i, , drop = FALSE], f = comp$f_orig))[1, ])
    }
  } else {
    stop("Invalid value of trans.")
  }
  res_score <- drop(res_score)
  return(res_score)
}
# ------------------------------------------------------------------------------


# Compute Fisher Information ---------------------------------------------------
#' @title Compute Fisher Information
#'
#' @description
#' A function computing Fisher information, its inverse, or its inverse square root for a given distribution.
#'
#' @inheritParams distr_density
#'
#' @param trans An optional transformation of the Fisher information. The supported transformations are the inverse of the Fisher information (\code{trans = "inv"}) and the inverse square root of the Fisher information (\code{trans = "inv_sqrt"}).
#'
#' @return The (transformed) Fisher information.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}}
#'
#' @examples
#' # Fisher information for the negative binomial distribution
#' distr_fisher(f = c(0.5, 1.2), distr = "negbin")
#'
#' # Fisher information for the multivariate normal distribution
#' distr_fisher(f = c(0, 0, 1, 1, 0.5), distr = "mnorm")
#'
#' @export
distr_fisher <- function(f, distr, param = NULL, par_link = NULL, trans = NULL) {
  comp <- list()
  comp$distr <- check_my_distr(distr = distr)
  comp$param <- check_my_param(param = param, distr = comp$distr)
  info_distr <- info_distribution(distr = comp$distr, param = comp$param)
  comp$f <- check_my_f(f = f)
  comp$t <- check_my_t(y = comp$f)
  comp$n <- check_my_n(f = comp$f, distr = comp$distr, param = comp$param, dim = info_distr$dim)
  info_par <- info_parameters(distr = comp$distr, param = comp$param, n = comp$n)
  comp$par_link <- check_my_par_link(par_link = par_link, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = comp$par_link)
  res_fisher <- array(NA_real_, dim = c(comp$t, info_par$par_num, info_par$par_num))
  if (is.null(trans)) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      res_fisher[i, , ] <- t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv
    }
  } else if (is.vector(trans) && !is.list(trans) && length(trans) == 1L && trans == "inv" && info_distr$orthog == TRUE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      res_fisher[i, , ] <- matrix_diag_inv(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv)
    }
  } else if (is.vector(trans) && !is.list(trans) && length(trans) == 1L && trans == "inv" && info_distr$orthog == FALSE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      res_fisher[i, , ] <- matrix_inv(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv)
    }
  } else if (is.vector(trans) && !is.list(trans) && length(trans) == 1L && trans == "inv_sqrt" && info_distr$orthog == TRUE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      res_fisher[i, , ] <- matrix_diag_inv_sqrt(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv)
    }
  } else if (is.vector(trans) && !is.list(trans) && length(trans) == 1L && trans == "inv_sqrt" && info_distr$orthog == FALSE) {
    for (i in 1:comp$t) {
      comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f[i, , drop = FALSE], par_trans = info_par$par_trans)
      comp$jacob_inv <- reparam_link_jacob_inv(f_orig = comp$f_orig, par_trans = info_par$par_trans)
      res_fisher[i, , ] <- matrix_inv_sqrt(t(comp$jacob_inv) %*% do.call(paste("distr", comp$distr, comp$param, "fisher", sep = "_"), args = list(f = comp$f_orig))[1, , ] %*% comp$jacob_inv)
    }
  } else {
    stop("Invalid value of trans.")
  }
  res_fisher <- drop(res_fisher)
  return(res_fisher)
}
# ------------------------------------------------------------------------------


# Generate Random Observations -------------------------------------------------
#' @title Generate Random Observations
#'
#' @description
#' A function generating random observations from a given distribution.
#'
#' @inheritParams distr_density
#'
#' @param t A number of generated observations.
#' @param f A numeric vector of parameters. The same parameters are used for each generated observation.
#'
#' @return The generated observations.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}}
#'
#' @examples
#' # Random observations from the negative binomial distribution
#' distr_random(t = 10, f = c(0.5, 1.2), distr = "negbin")
#'
#' # Random observations from the multivariate normal distribution
#' distr_random(t = 10, f = c(0, 0, 1, 1, 0.5), distr = "mnorm")
#'
#' @export
distr_random <- function(t, f, distr, param = NULL, par_link = NULL) {
  comp <- list()
  comp$distr <- check_my_distr(distr = distr)
  comp$param <- check_my_param(param = param, distr = comp$distr)
  info_distr <- info_distribution(distr = comp$distr, param = comp$param)
  comp$f <- check_my_f(f = f, t = 1L)
  comp$t <- check_my_t(t = t)
  comp$n <- check_my_n(f = comp$f, distr = comp$distr, param = comp$param, dim = info_distr$dim)
  info_par <- info_parameters(distr = comp$distr, param = comp$param, n = comp$n)
  comp$par_link <- check_my_par_link(par_link = par_link, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  info_par <- info_linked_parameters(info_par = info_par, par_link = comp$par_link)
  comp$f_orig <- reparam_tilde_to_orig(f_tilde = comp$f, par_trans = info_par$par_trans)
  res_random <- do.call(paste("distr", comp$distr, comp$param, "random", sep = "_"), args = list(t = comp$t, f = comp$f_orig))
  res_random <- drop(res_random)
  return(res_random)
}
# ------------------------------------------------------------------------------


