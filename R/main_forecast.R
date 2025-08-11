
# FORECASTING OF GAS MODEL FUNCTION


# Forecast GAS Model -----------------------------------------------------------
#' @title Forecast GAS Model
#'
#' @description
#' A function for forecasting of generalized autoregressive score (GAS) models of Creal et al. (2013) and Harvey (2013).
#' Method \code{"mean_path"} filters time-varying parameters based on zero score and then generates mean of time series.
#' Method \code{"simulated_paths"} repeatedly simulates time series, simultaneously filters time-varying parameters, and then estimates mean, standard deviation, and quantiles (see Blasques et al., 2016).
#' Instead of supplying arguments about the model, the function can be applied to the \code{gas} object obtained by the \code{\link[gasmodel:gas]{gas()}} function.
#'
#' @param gas_object An optional GAS estimate, i.e. a list of S3 class \code{gas} returned by function \code{\link[gasmodel:gas]{gas()}}.
#' @param method  A method used for forecasting. Supported methods are \code{"mean_path"} and \code{"simulated_paths"}.
#' @param t_ahead A number of observations to forecast.
#' @param x_ahead Out-of-sample exogenous variables. For a single variable common for all time-varying parameters, a numeric vector. For multiple variables common for all time-varying parameters, a numeric matrix with observations in rows. For individual variables for each time-varying parameter, a list of numeric vectors or matrices in the above form. The number of observation must be equal to \code{t_ahead}.
#' @param rep_ahead A number of simulation repetitions for \code{method = "simulated_paths"}.
#' @param quant A numeric vector of probabilities determining quantiles for \code{method = "simulated_paths"}.
#' @param y,x,distr,param,scaling,regress,p,q,par_static,par_link,par_init,coef_est When \code{gas_object} is not supplied, the estimated model can be specified using these individual arguments. See the arguments and value of the \code{\link[gasmodel:gas]{gas()}} function for more details.
#'
#' @return A \code{list} of S3 class \code{gas_forecast} with components:
#' \item{data$y}{The time series.}
#' \item{data$x}{The exogenous variables.}
#' \item{data$x_ahead}{The out-of-sample exogenous variables.}
#' \item{model$distr}{The conditional distribution.}
#' \item{model$param}{The parametrization of the conditional distribution.}
#' \item{model$scaling}{The scaling function.}
#' \item{model$regress}{The specification of the regression and dynamic equation.}
#' \item{model$t}{The length of the time series.}
#' \item{model$t_ahead}{The length of the out-of-sample time series.}
#' \item{model$n}{The dimension of the model.}
#' \item{model$m}{The number of exogenous variables.}
#' \item{model$p}{The score order.}
#' \item{model$q}{The autoregressive order.}
#' \item{model$par_static}{The static parameters.}
#' \item{model$par_link}{The parameters with the logarithmic/logistic links.}
#' \item{model$par_init}{The initial values of the time-varying parameters.}
#' \item{model$coef_est}{The estimated coefficients.}
#' \item{forecast$method}{The method used for forecasting.}
#' \item{forecast$y_ahead_mean}{The mean of the forecasted time series.}
#' \item{forecast$y_ahead_sd}{The standard deviation of the forecasted time series. Only for \code{method = "simulated_paths"}.}
#' \item{forecast$y_ahead_quant}{The quantiles of the forecasted time series. Only for \code{method = "simulated_paths"}.}
#' \item{forecast$par_tv_ahead_mean}{The mean of the forecasted time-varying parameters.}
#' \item{forecast$par_tv_ahead_sd}{The standard deviation of the forecasted time-varying parameters. Only for \code{method = "simulated_paths"}.}
#' \item{forecast$par_tv_ahead_quant}{The quantiles of the forecasted time-varying parameters. Only for \code{method = "simulated_paths"}.}
#' \item{forecast$score_tv_ahead_mean}{The mean of the forecasted scores.}
#' \item{forecast$score_tv_ahead_sd}{The standard deviation of the forecasted scores. Only for \code{method = "simulated_paths"}.}
#' \item{forecast$score_tv_ahead_quant}{The quantiles of the forecasted scores. Only for \code{method = "simulated_paths"}.}
#'
#' @note
#' Supported generic functions for S3 class \code{gas_forecast} include \code{\link[base:summary]{summary()}} ans \code{\link[base:plot]{plot()}}.
#'
#' @references
#' Blasques, F., Koopman, S. J., Łasak, K., and Lucas, A. (2016). In-Sample Confidence Bands and Out-of-Sample Forecast Bands for Time-Varying Parameters in Observation-Driven Models. \emph{International Journal of Forecasting}, \strong{32}(3), 875–887. \doi{10.1016/j.ijforecast.2015.11.018}.
#'
#' Creal, D., Koopman, S. J., and Lucas, A. (2013). Generalized Autoregressive Score Models with Applications. \emph{Journal of Applied Econometrics}, \strong{28}(5), 777–795. \doi{10.1002/jae.1279}.
#'
#' Harvey, A. C. (2013). \emph{Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series}. Cambridge University Press. \doi{10.1017/cbo9781139540933}.
#'
#' @seealso
#' \code{\link[gasmodel:gas]{gas()}}
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
#' # Forecast the model by the "mean_paths" method
#' x_ahead <- cbind(kronecker(matrix(1, 53, 1), diag(7)), 1)[3:367, -1]
#' fcst_negbin <- gas_forecast(est_negbin, t_ahead = 365, x_ahead = x_ahead)
#' fcst_negbin
#'
#' # Plot the forecasted expected value
#' plot(fcst_negbin)}
#'
#' @export
gas_forecast <- function(gas_object = NULL, method = "mean_path", t_ahead = 1L, x_ahead = NULL, rep_ahead = 1000L, quant = c(0.025, 0.975), y = NULL, x = NULL, distr = NULL, param = NULL, scaling = "unit", regress = "joint", p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, coef_est = NULL) {
  if (!is.null(gas_object) && inherits(gas_object, "gas")) {
    gas_forecast(gas_object = NULL, method = method, t_ahead = t_ahead, x_ahead = x_ahead, rep_ahead = rep_ahead, quant = quant, y = gas_object$data$y, x = gas_object$data$x, distr = gas_object$model$distr, param = gas_object$model$param, scaling = gas_object$model$scaling, regress = gas_object$model$regress, p = gas_object$model$p, q = gas_object$model$q, par_static = gas_object$model$par_static, par_link = gas_object$model$par_link, par_init = gas_object$model$par_init, coef_est = gas_object$fit$coef_est)
  } else if (!is.null(gas_object)) {
    stop("Unsupported class of gas_object.")
  } else {
    # Load auxiliary variables:
    load <- load_forecast(method = method, t_ahead = t_ahead, x_ahead = x_ahead, rep_ahead = rep_ahead, quant = quant, y = y, x = x, distr = distr, param = param, scaling = scaling, regress = regress, p = p, q = q, par_static = par_static, par_link = par_link, par_init = par_init, coef_est = coef_est)
    data <- load$data
    model <- load$model
    fun <- load$fun
    info_distr <- load$info_distr
    info_par <- load$info_par
    info_coef <- load$info_coef
    comp <- load$comp
    forecast <- load$forecast
    # Forecast using mean path method:
    if (forecast$method == "mean_path") {
      # Compute path for static model:
      if (all(model$p + model$q == 0L)) {
        comp$par_init <- model$par_init
        if (any(is.na(comp$par_init))) {
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
        for (j in comp$idx_ok_regular) {
          comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
        }
        comp$score_tv[comp$idx_ok_ahead, ] <- 0
      # Compute path for dynamic joint model:
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
        for (j in comp$idx_ok_regular) {
          for (i in 1:info_par$par_num) {
            cur_e[i] <- sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
          }
          comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
          comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
        }
        comp$score_tv[comp$idx_ok_ahead, ] <- 0
        for (j in comp$idx_ok_ahead) {
          for (i in 1:info_par$par_num) {
            cur_e[i] <- sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
          }
          comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
        }
      # Compute path for dynamic separate model:
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
        for (j in comp$idx_ok_regular) {
          for (i in 1:info_par$par_num) {
            cur_e[i] <- sum(comp$err_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
          }
          comp$err_tv[j, ] <- cur_e
          comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
          comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
        }
        comp$score_tv[comp$idx_ok_ahead, ] <- 0
        for (j in comp$idx_ok_ahead) {
          for (i in 1:info_par$par_num) {
            cur_e[i] <- sum(comp$err_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
          }
          comp$err_tv[j, ] <- cur_e
        }
        comp$par_tv[comp$idx_ok_ahead, ] <- comp$par_tv[comp$idx_ok_ahead, ] + comp$err_tv[comp$idx_ok_ahead, ]
      }
      # Format results:
      comp$y[comp$idx_ok_ahead, ] <- fun$mean(f = comp$par_tv[comp$idx_ok_ahead, , drop = FALSE])
      info_data <- info_data(y = data$y, x = data$x)
      data$y <- name_matrix(data$y, info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
      data$x <- name_list_of_matrices(data$x, info_par$par_names, info_data$index_time_list, info_data$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
      info_data_ahead <- info_data(y = matrix(nrow = model$t_ahead, ncol = model$n), x = data$x_ahead, skip_t = model$t)
      data$x_ahead <- name_list_of_matrices(data$x_ahead, info_par$par_names, info_data_ahead$index_time_list, info_data_ahead$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
      forecast$y_ahead_mean <- name_matrix(comp$y[(comp$pre_num + model$t + 1L):comp$full_num, , drop = FALSE], info_data_ahead$index_time, info_data_ahead$index_series, drop = c(FALSE, TRUE))
      forecast$par_tv_ahead_mean <- name_matrix(comp$par_tv[(comp$pre_num + model$t + 1L):comp$full_num, , drop = FALSE], info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      forecast$score_tv_ahead_mean <- name_matrix(comp$score_tv[(comp$pre_num + model$t + 1L):comp$full_num, , drop = FALSE], info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
    # Forecast using simulated path method:
    } else if (forecast$method == "simulated_paths") {
      comp$rep_ahead <- check_generic_positive_integer_scalar(arg = rep_ahead, arg_name = "rep_ahead")
      comp$quant <- check_generic_probability_vector(arg = quant, arg_name = "quant")
      comp$y_ahead_path <- array(NA_real_, dim = c(comp$rep_ahead, model$t_ahead, model$n))
      comp$par_tv_ahead_path <- array(NA_real_, dim = c(comp$rep_ahead, model$t_ahead, info_par$par_num))
      comp$score_tv_ahead_path <- array(NA_real_, dim = c(comp$rep_ahead, model$t_ahead, info_par$par_num))
      # Compute path for static model:
      if (all(model$p + model$q == 0L)) {
        comp$par_init <- model$par_init
        if (any(is.na(comp$par_init))) {
          comp$par_unc <- sapply(1:info_par$par_num, function(i) { (comp$omega_vector[i] + comp$average_x[[i]] %*% comp$beta_list[[i]]) })
          comp$par_init[is.na(comp$par_init)] <- comp$par_unc[is.na(comp$par_init)]
        }
        if (length(comp$idx_na) > 0L) {
          comp$par_tv[comp$idx_na, ] <- matrix(comp$par_init, nrow = length(comp$idx_na), ncol = info_par$par_num, byrow = TRUE)
          comp$score_tv[comp$idx_na, ] <- 0
        }
        comp$par_tv[comp$idx_ok, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok), ncol = info_par$par_num, byrow = TRUE)
        if (any(model$m > 0L)) {
          comp$par_tv[comp$idx_ok, ] <- comp$par_tv[comp$idx_ok, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_regular, , drop = FALSE] %*% comp$beta_list[[i]] })
        }
        for (j in comp$idx_ok_regular) {
          comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
        }
        for (a in 1:comp$rep_ahead) {
          for (j in comp$idx_ok_ahead) {
            comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
            comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
          }
          comp$y_ahead_path[a, , ] <- comp$y[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$par_tv_ahead_path[a, , ] <- comp$par_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$score_tv_ahead_path[a, , ] <- comp$score_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
        }
      # Compute path for dynamic joint model:
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
        comp$par_tv[comp$idx_ok_regular, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok_regular), ncol = info_par$par_num, byrow = TRUE)
        if (any(model$m > 0L)) {
          comp$par_tv[comp$idx_ok_regular, ] <- comp$par_tv[comp$idx_ok_regular, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_regular, , drop = FALSE] %*% comp$beta_list[[i]] })
        }
        cur_e <- rep(NA_real_, info_par$par_num)
        for (j in comp$idx_ok_regular) {
          for (i in 1:info_par$par_num) {
            cur_e[i] <- sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
          }
          comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
          comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
        }
        for (a in 1:comp$rep_ahead) {
          comp$par_tv[comp$idx_ok_ahead, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok_ahead), ncol = info_par$par_num, byrow = TRUE)
          if (any(model$m > 0L)) {
            comp$par_tv[comp$idx_ok_ahead, ] <- comp$par_tv[comp$idx_ok_ahead, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_ahead, , drop = FALSE] %*% comp$beta_list[[i]] })
          }
          cur_e <- rep(NA_real_, info_par$par_num)
          for (j in comp$idx_ok_ahead) {
            for (i in 1:info_par$par_num) {
              cur_e[i] <- sum(comp$par_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
            }
            comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
            comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
            comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
          }
          comp$y_ahead_path[a, , ] <- comp$y[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$par_tv_ahead_path[a, , ] <- comp$par_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$score_tv_ahead_path[a, , ] <- comp$score_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
        }
      # Compute path for dynamic separate model:
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
        comp$par_tv[comp$idx_ok_regular, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok_regular), ncol = info_par$par_num, byrow = TRUE)
        if (any(model$m > 0L)) {
          comp$par_tv[comp$idx_ok_regular, ] <- comp$par_tv[comp$idx_ok_regular, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_regular, , drop = FALSE] %*% comp$beta_list[[i]] })
        }
        cur_e <- rep(NA_real_, info_par$par_num)
        for (j in comp$idx_ok_regular) {
          for (i in 1:info_par$par_num) {
            cur_e[i] <- sum(comp$err_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
          }
          comp$err_tv[j, ] <- cur_e
          comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
          comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
        }
        for (a in 1:comp$rep_ahead) {
          comp$par_tv[comp$idx_ok_ahead, ] <- matrix(comp$omega_vector, nrow = length(comp$idx_ok_ahead), ncol = info_par$par_num, byrow = TRUE)
          if (any(model$m > 0L)) {
            comp$par_tv[comp$idx_ok_ahead, ] <- comp$par_tv[comp$idx_ok_ahead, ] + sapply(1L:info_par$par_num, function(i) { comp$x[[i]][comp$idx_ok_ahead, , drop = FALSE] %*% comp$beta_list[[i]] })
          }
          cur_e <- rep(NA_real_, info_par$par_num)
          for (j in comp$idx_ok_ahead) {
            for (i in 1:info_par$par_num) {
              cur_e[i] <- sum(comp$err_tv[j - seq_along(comp$phi_list[[i]]), i] * comp$phi_list[[i]]) + sum(comp$score_tv[j - seq_along(comp$alpha_list[[i]]), i] * comp$alpha_list[[i]])
            }
            comp$err_tv[j, ] <- cur_e
            comp$par_tv[j, ] <- comp$par_tv[j, ] + cur_e
            comp$y[j, ] <- fun$random(t = 1L, f = comp$par_tv[j, , drop = FALSE])
            comp$score_tv[j, ] <- fun$score(y = comp$y[j, , drop = FALSE], f = comp$par_tv[j, , drop = FALSE])
          }
          comp$y_ahead_path[a, , ] <- comp$y[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$par_tv_ahead_path[a, , ] <- comp$par_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
          comp$score_tv_ahead_path[a, , ] <- comp$score_tv[(comp$pre_num + model$t + 1L):comp$full_num, ]
        }
      }
      # Format results:
      info_data <- info_data(y = data$y, x = data$x)
      data$y <- name_matrix(data$y, info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
      data$x <- name_list_of_matrices(data$x, info_par$par_names, info_data$index_time_list, info_data$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
      info_data_ahead <- info_data(y = matrix(nrow = model$t_ahead, ncol = model$n), x = data$x_ahead, skip_t = model$t)
      data$x_ahead <- name_list_of_matrices(data$x_ahead, info_par$par_names, info_data_ahead$index_time_list, info_data_ahead$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
      forecast$y_ahead_mean <- name_matrix(colMeans(comp$y_ahead_path, na.rm = TRUE), info_data_ahead$index_time, info_data_ahead$index_series, drop = c(FALSE, TRUE))
      forecast$y_ahead_sd <- name_matrix(apply(comp$y_ahead_path, 2:3, stats::sd, na.rm = TRUE), info_data_ahead$index_time, info_data_ahead$index_series, drop = c(FALSE, TRUE))
      forecast$y_ahead_quant <- name_3d_array(aperm(array(apply(comp$y_ahead_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t_ahead, model$n)), c(2, 3, 1)), info_data_ahead$index_time, info_data_ahead$index_series, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
      forecast$par_tv_ahead_mean <- name_matrix(colMeans(comp$par_tv_ahead_path, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      forecast$par_tv_ahead_sd <- name_matrix(apply(comp$par_tv_ahead_path, 2:3, stats::sd, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      forecast$par_tv_ahead_quant <- name_3d_array(aperm(array(apply(comp$par_tv_ahead_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t_ahead, info_par$par_num)), c(2, 3, 1)), info_data_ahead$index_time, info_par$par_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
      forecast$score_tv_ahead_mean <- name_matrix(colMeans(comp$score_tv_ahead_path, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      forecast$score_tv_ahead_sd <- name_matrix(apply(comp$score_tv_ahead_path, 2:3, stats::sd, na.rm = TRUE), info_data_ahead$index_time, info_par$par_names, drop = c(FALSE, TRUE))
      forecast$score_tv_ahead_quant <- name_3d_array(aperm(array(apply(comp$score_tv_ahead_path, 2:3, stats::quantile, probs = comp$quant, na.rm = TRUE), dim = c(length(comp$quant), model$t_ahead, info_par$par_num)), c(2, 3, 1)), info_data_ahead$index_time, info_par$par_names, paste0(comp$quant * 100, "%"), drop = c(FALSE, TRUE, TRUE), zero = c(FALSE, FALSE, TRUE))
    }
    # Report results
    report <- list(data = data, model = model, forecast = forecast)
    class(report) <- "gas_forecast"
    return(report)
  }
}
# ------------------------------------------------------------------------------


# Print Forecasts --------------------------------------------------------------
#' @export
print.gas_forecast <- function(x, ...) {
  info_title <- info_title(distr = x$model$distr, param = x$model$param, scaling = x$model$scaling)
  cat("GAS Model:", info_title$title, "\n")
  cat("\n")
  cat("Method:", switch(x$forecast$method, "mean_path" = "Mean Path", "simulated_paths" = "Simulated Paths"), "\n")
  cat("\n")
  cat("Forecasts: \n")
  print(x$forecast$y_ahead_mean)
  invisible(x)
}
# ------------------------------------------------------------------------------


# Summarize Forecasts ----------------------------------------------------------
#' @export
summary.gas_forecast <- function(object, ...) {
  print(object)
  cat("\n")
  cat("Time-Varying Parameters:", "\n")
  print(object$forecast$par_tv_ahead_mean)
  invisible(object)
}
# ------------------------------------------------------------------------------


# Plot Forecasted Time Series ---------------------------------------------------
#' @importFrom dplyr %>%
#' @importFrom ggplot2 .data
#' @export
plot.gas_forecast <- function(x, which = NULL, ...) {
  y <- x$data$y
  y_fcst <- x$forecast$y_ahead_mean
  if (is.vector(y)) {
    y_full <- c(y, y_fcst)
    ts_index <- 1:length(y_full)
    ts_divide <- length(y) + 0.5
    gg_data <- dplyr::tibble(index = ts_index, value = y_full)
    gg_fig <- ggplot2::ggplot(gg_data, mapping = ggplot2::aes(.data$index, .data$value)) +
      ggplot2::geom_vline(xintercept = ts_divide, linetype = "dotted") +
      ggplot2::geom_line(color = "#800000") +
      ggplot2::geom_point(color = "#800000") +
      ggplot2::labs(title = "Forecasted Time Series", x = "Time Index", y = "Observation Value")
    gg_list <- list(gg_fig)
  } else {
    y_full <- rbind(y, y_fcst)
    ser_names <- colnames(y_full)
    ts_index <- 1:nrow(y_full)
    ts_divide <- nrow(y) + 0.5
    gg_col <- c(rep("#CCCCCC", times = ncol(y_full) - 1), "#800000")
    gg_list <- list()
    for (i in 1:length(ser_names)) {
      gg_levels <- c(ser_names[-i], ser_names[i])
      gg_data <-  y_full %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(index = ts_index) %>%
        tidyr::pivot_longer(cols = -dplyr::last_col(), names_to = "ser", values_to = "value") %>%
        dplyr::mutate(ser = factor(.data$ser, levels = gg_levels, ordered = TRUE)) %>%
        dplyr::arrange(.data$ser)
      gg_fig <- ggplot2::ggplot(gg_data, mapping = ggplot2::aes(.data$index, .data$value, group = .data$ser, color = .data$ser)) +
        ggplot2::geom_vline(xintercept = ts_divide, linetype = "dotted") +
        ggplot2::geom_line(show.legend = FALSE) +
        ggplot2::geom_point(show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = gg_col) +
        ggplot2::labs(title = paste("Forecasted Time Series", ser_names[i]), x = "Time Index", y = "Observation Value")
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


