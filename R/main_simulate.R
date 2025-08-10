
# SIMULATION OF GAS MODEL FUNCTION


# Simulate GAS Model -----------------------------------------------------------
#' @title Simulate GAS Model
#'
#' @description
#' A function for simulation of generalized autoregressive score (GAS) models of Creal et al. (2013) and Harvey (2013).
#' Instead of supplying arguments about the model, the function can be applied to the \code{gas} object obtained by the \code{\link[gasmodel:gas]{gas()}} function.
#'
#' @param gas_object An optional GAS estimate, i.e. a list of S3 class \code{gas} returned by function \code{\link[gasmodel:gas]{gas()}}.
#' @param t_sim A number of observations to simulate.
#' @param x_sim Exogenous variables used for simulations. For a single variable common for all time-varying parameters, a numeric vector. For multiple variables common for all time-varying parameters, a numeric matrix with observations in rows. For individual variables for each time-varying parameter, a list of numeric vectors or matrices in the above form. The number of observation must be equal to \code{t_sim}.
#' @param distr,param,scaling,regress,n,p,q,par_static,par_link,par_init,coef_est When \code{gas_object} is not supplied, the estimated model can be specified using these individual arguments. See the arguments and value of the \code{\link[gasmodel:gas]{gas()}} function for more details.
#'
#' @return A \code{list} of S3 class \code{gas_simulate} with components:
#' \item{data$x_sim}{The exogenous variables used in simulation.}
#' \item{model$distr}{The conditional distribution.}
#' \item{model$param}{The parametrization of the conditional distribution.}
#' \item{model$scaling}{The scaling function.}
#' \item{model$regress}{The specification of the regression and dynamic equation.}
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
#' @note
#' Supported generic functions for S3 class \code{gas_simulate} include \code{\link[base:summary]{summary()}} ans \code{\link[base:plot]{plot()}}.
#'
#' @references
#' Creal, D., Koopman, S. J., and Lucas, A. (2013). Generalized Autoregressive Score Models with Applications. \emph{Journal of Applied Econometrics}, \strong{28}(5), 777–795. \doi{10.1002/jae.1279}.
#'
#' Harvey, A. C. (2013). \emph{Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series}. Cambridge University Press. \doi{10.1017/cbo9781139540933}.
#'
#' @seealso \code{\link[gasmodel:gas]{gas()}}
#'
#' @examples
#' \donttest{# Simulate GAS model based on the negative binomial distribution
#' sim_negbin <- gas_simulate(t_sim = 50, distr = "negbin", reg = "sep",
#'   coef_est = c(2.60, 0.02, 0.95, 0.03))
#' sim_negbin
#'
#' # Plot the simulated time series
#' plot(sim_negbin)}
#'
#' @export
gas_simulate <- function(gas_object = NULL, t_sim = 1L, x_sim = NULL, distr = NULL, param = NULL, scaling = "unit", regress = "joint", n = NULL, p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, coef_est = NULL) {
  if (!is.null(gas_object) && inherits(gas_object, "gas")) {
    gas_simulate(gas_object = NULL, t_sim = t_sim, x_sim = x_sim, distr = gas_object$model$distr, param = gas_object$model$param, scaling = gas_object$model$scaling, regress = gas_object$model$regress, n = gas_object$model$n, p = gas_object$model$p, q = gas_object$model$q, par_static = gas_object$model$par_static, par_link = gas_object$model$par_link, par_init = gas_object$model$par_init, coef_est = gas_object$fit$coef_est)
  } else if (!is.null(gas_object)) {
    stop("Unsupported class of gas_object.")
  } else {
    load <- load_simulate(t_sim = t_sim, x_sim = x_sim, distr = distr, param = param, scaling = scaling, regress = regress, n = n, p = p, q = q, par_static = par_static, par_link = par_link, par_init = par_init, coef_est = coef_est)
    data <- load$data
    model <- load$model
    fun <- load$fun
    info_distr <- load$info_distr
    info_par <- load$info_par
    info_coef <- load$info_coef
    comp <- load$comp
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
    } else if (model$regress == "joint") {
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
    } else if (model$regress == "sep") {
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


# Summarize Simulations --------------------------------------------------------
#' @export
summary.gas_simulate <- function(object, ...) {
  print(object)
  cat("\n")
  cat("Time-Varying Parameters:", "\n")
  print(object$simulation$par_tv_sim)
  invisible(object)
}
# ------------------------------------------------------------------------------


# Plot Simulated Time Series ---------------------------------------------------
#' @importFrom dplyr %>%
#' @importFrom ggplot2 .data
#' @export
plot.gas_simulate <- function(x, which = NULL, ...) {
  y_sim <- x$simulation$y_sim
  if (is.vector(y_sim)) {
    ts_index <- 1:length(y_sim)
    gg_data <- dplyr::tibble(index = ts_index, value = y_sim)
    gg_fig <- ggplot2::ggplot(gg_data, mapping = ggplot2::aes(.data$index, .data$value)) +
      ggplot2::geom_line(color = "#800000") +
      ggplot2::geom_point(color = "#800000") +
      ggplot2::labs(title = "Simulated Time Series", x = "Time Index", y = "Observation Value")
    gg_list <- list(gg_fig)
  } else {
    ser_names <- colnames(y_sim)
    ts_index <- 1:nrow(y_sim)
    gg_col <- c(rep("#CCCCCC", times = ncol(y_sim) - 1), "#800000")
    gg_list <- list()
    for (i in 1:length(ser_names)) {
      gg_levels <- c(ser_names[-i], ser_names[i])
      gg_data <-  y_sim %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(index = ts_index) %>%
        tidyr::pivot_longer(cols = -dplyr::last_col(), names_to = "ser", values_to = "value") %>%
        dplyr::mutate(ser = factor(.data$ser, levels = gg_levels, ordered = TRUE)) %>%
        dplyr::arrange(.data$ser)
      gg_fig <- ggplot2::ggplot(gg_data, mapping = ggplot2::aes(.data$index, .data$value, group = .data$ser, color = .data$ser)) +
        ggplot2::geom_line(show.legend = FALSE) +
        ggplot2::geom_point(show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = gg_col) +
        ggplot2::labs(title = paste("Simulated Time Series", ser_names[i]), x = "Time Index", y = "Observation Value")
      gg_list <- append(gg_list, list(gg_fig))
    }
  }
  gg_which <- 1:length(gg_list)
  if (!is.null(which)) {
    gg_which <- gg_which[gg_which %in% which]
  }
  if (length(gg_which) == 1) {
    be_silent(print(gg_list[[gg_which[1]]]))
  } else if (length(gg_which) > 1) {
    be_silent(print(gg_list[[gg_which[1]]]))
    old_par <- grDevices::devAskNewPage(ask = TRUE)
    for (i in 2:length(gg_which)) {
      be_silent(print(gg_list[[gg_which[i]]]))
    }
    on.exit(grDevices::devAskNewPage(ask = old_par))
  }
  invisible(gg_list)
}
# ------------------------------------------------------------------------------


