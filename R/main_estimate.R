
# ESTIMATION OF GAS MODEL FUNCTION


# Estimate GAS Model -----------------------------------------------------------
#' @title Estimate GAS Model
#'
#' @description
#' A versatile function for estimation of generalized autoregressive score (GAS) models of Creal et al. (2013) and Harvey (2013).
#' Model specification allows for various conditional distributions, different parametrizations, exogenous variables, higher score and autoregressive orders, custom and unconditional initial values of time-varying parameters, fixed and bounded values of coefficients, and NA values.
#' Model estimation is performed by the maximum likelihood method and the Hessian matrix.
#' The function can be supplied with any optimization and Hessian functions.
#'
#' @param y A time series. For univariate time series, a numeric vector or a matrix with a single column. For multivariate times series, a numeric matrix with observations in rows.
#' @param x Optional exogenous variables. For a single variable common for all time-varying parameters, a numeric vector. For multiple variables common for all time-varying parameters, a numeric matrix with observations in rows. For individual variables for each time-varying parameter, a list of numeric vectors or matrices in the above form. The number of observation must be equal to the number of observations of \code{y}.
#' @param distr A conditional distribution. See \code{\link[gasmodel:distr]{distr()}} for available distributions.
#' @param param A parametrization of the conditional distribution. If \code{NULL}, default parametrization is used. See \code{\link[gasmodel:distr]{distr()}} for available parametrizations.
#' @param scaling A scaling function for the score. The supported scaling functions are the unit scaling (\code{scaling = "unit"}), the inverse of the Fisher information matrix scaling (\code{scaling = "fisher_inv"}), and the inverse square root of the Fisher information matrix scaling (\code{scaling = "fisher_inv_sqrt"}). The latter two scalings use the Fisher information for the time-varying parameters only. For the full Fisher information matrix for both time-varying and static parameters, there are the \code{"full_fisher_inv"} and \code{"full_fisher_inv_sqrt"} scalings. For the individual Fisher information for each parameter, there are the \code{"diag_fisher_inv"} and \code{"diag_fisher_inv_sqrt"} scalings. Note that when the parametrization is orthogonal (see \code{\link[gasmodel:distr]{distr()}}), there are no differences between these scaling variants.
#' @param regress A specification of the regression and dynamic equation with regard to exogenous variables. The supported specifications are exogenous variables and dynamics within the same equation (\code{regress = "joint"}) and separate equations for exogenous variables and dynamics in the fashion of regression models with dynamic errors (\code{regress = "sep"}). In a stationary model without exogenous variables, the two specifications are equivalent, although with differently parametrized intercept.
#' @param p A score order. For order common for all parameters, a numeric vector of length 1. For individual order for each parameter, a numeric vector of length equal to the number of parameters. Defaults to \code{1L}.
#' @param q An autoregressive order. For order common for all parameters, a numeric vector of length 1. For individual order for each parameter, a numeric vector of length equal to the number of parameters. Defaults to \code{1L}.
#' @param par_static An optional logical vector indicating static parameters. Overrides \code{x}, \code{p}, and \code{q}.
#' @param par_link An optional logical vector indicating whether the logarithmic/logistic link should be applied to restricted parameters in order to obtain unrestricted values. Defaults to applying the logarithmic/logistic link for time-varying parameters and keeping the original link for constant parameters.
#' @param par_init An optional numeric vector of initial values of time-varying parameters. For \code{NA} values or when \code{NULL}, set initial values to unconditional values of time-varying parameters. For example, in the case of GAS(1,1) model with \code{regress = "joint"}, to \code{omega / (1 - phi1)}. Not to be confused with starting values for the optimization \code{coef_start}.
#' @param lik_skip A numeric value specifying the number of skipped observations at the beginning of the time series or after \code{NA} values in the likelihood computation. Defaults to \code{0L}, i.e. the full likelihood. If \code{NULL}, it is selected as \code{max(p,q)}, i.e. the conditional likelihood.
#' @param coef_fix_value An optional numeric vector of values to which coefficients are to be fixed. \code{NA} values represent coefficients to be estimated.
#' @param coef_fix_other An optional square numeric matrix of multiples of the estimated coefficients, which are to be added to the fixed coefficients. This allows the fixed coefficients to be linear combinations of the estimated coefficients. A coefficient given by row is fixed on coefficient given by column. By this logic, all rows corresponding to the estimated coefficients should contain only \code{NA} values. Furthermore, all columns corresponding to the fixed coefficients should also contain only \code{NA} values.
#' @param coef_fix_special An optional character vector of predefined structures of \code{coef_fix_value} and \code{coef_fix_other}. Useful mainly for multidimensional models. Value \code{"panel_structure"} forces all regression, autoregression, and score coefficients to be the same for all time-varying parameters within their group. Value \code{"zero_sum_intercept"} forces all constant parameters to sum up to zero within their group. Value \code{"random_walk"} forces all autoregressive coefficients to be equal to one (should be used with caution due to nonstationarity; \code{par_init} must be specified). Multiple predefined structures can be used together. Also can be used in combination with custom \code{coef_fix_value} and \code{coef_fix_other}.
#' @param coef_bound_lower An optional numeric vector of lower bounds on coefficients.
#' @param coef_bound_upper An optional numeric vector of upper bounds on coefficients.
#' @param coef_start An optional numeric vector of starting values for the optimization. If not supplied, starting values are selected from a small grid of values.
#' @param optim_function An optimization function. For suitable wrappers of common R optimization functions, see \code{\link[gasmodel:wrappers_optim]{wrappers_optim}}. Can be set to \code{NULL} if the optimal solution should not be computed, which can be useful if the goal is only to evaluate the fit for the coefficients specified in argument \code{coef_start}.
#' @param optim_arguments An optional list of arguments to be passed to the optimization function.
#' @param hessian_function A Hessian function. For suitable wrappers of common R Hessian functions, see \code{\link[gasmodel:wrappers_hessian]{wrappers_hessian}}. Can be set to \code{NULL} if the Hessian matrix should not be computed, which can speed up computations when asymptotic inference is not desired.
#' @param hessian_arguments An optional list of arguments to be passed to the Hessian function.
#' @param print_progress A logical value indicating whether to progressively print a detailed report on computation.
#'
#' @details
#' The generalized autoregressive score (GAS) models of Creal et al. (2013) and Harvey (2013), also known as dynamic conditional score (DCS) models or score-driven (SD) models, have established themselves as a useful modern framework for time series modeling.
#'
#' The GAS models are observation-driven models allowing for any underlying probability distribution \eqn{p(y_t|f_t)} with any time-varying parameters \eqn{f_t} for time series \eqn{y_t}.
#' They capture the dynamics of time-varying parameters using the autoregressive term and the lagged score, i.e. the gradient of the log-likelihood function.
#' Exogenous variables can also be included.
#' Specifically, time-varying parameters \eqn{f_t} follow the recursion
#' \deqn{f_t = \omega + \sum_{i=1}^M \beta_i x_{ti} + \sum_{j=1}^P \alpha_j S(f_{t-j}) \nabla(y_{t-j}, f_{t-j}) + \sum_{k=1}^Q \varphi_k f_{t-k},}{f_t = \omega + \sum_{i=1}^M \beta_i x_{ti} + \sum_{j=1}^P \alpha_j S(f_{t-j}) ∇(y_{t-j}, f_{t-j}) + \sum_{k=1}^Q \phi_k f_{t-k},}
#' where \eqn{\omega} is the intercept, \eqn{\beta_i} are the regression parameters, \eqn{\alpha_j} are the score parameters, \eqn{\varphi_k}{\phi_k} are the autoregressive parameters, \eqn{x_{ti}} are the exogenous variables, \eqn{S(f_t)} is a scaling function for the score, and \eqn{\nabla(y_t, f_t)}{∇(y_t, f_t)} is the score given by
#' \deqn{\nabla(y_t, f_t) = \frac{\partial \ln p(y_t | f_t)}{\partial f_t}.}{∇(y_t, f_t) = d ln p(y_t | f_t) / d f_t.}
#' In the case of a single time-varying parameter, \eqn{\omega}, \eqn{\beta_i}, \eqn{\alpha_j}, \eqn{\varphi_k}, \eqn{x_{ti}}, \eqn{S(f_t)}, and \eqn{\nabla(y_t, f_t)} are all scalar.
#' In the case of multiple time-varying parameters, \eqn{x_{ti}} are scalar, \eqn{\omega}, \eqn{\beta_i}, and \eqn{\nabla(y_{t - j}, f_{t - j})} are vectors, \eqn{\alpha_j} and \eqn{\varphi_k} are diagonal matrices, and \eqn{S(f_t)} is a square matrix.
#' Alternatively, a different model can be obtained by defining the recursion in the fashion of regression models with dynamic errors as
#' \deqn{f_t = \omega + \sum_{i=1}^M \beta_i x_{ti} + e_{t}, \quad e_t = \sum_{j=1}^P \alpha_j S(f_{t-j}) \nabla(y_{t-j}, f_{t-j}) + \sum_{k=1}^Q \varphi_k e_{t-k}.}{f_t = \omega + \sum_{i=1}^M \beta_i x_{ti} + e_t,   e_t = \sum_{j=1}^P \alpha_j S(f_{t-j}) ∇(y_{t-j}, f_{t-j}) + \sum_{k=1}^Q \phi_k e_{t-k}.}
#'
#' The GAS models can be straightforwardly estimated by the maximum likelihood method.
#' For the asymptotic theory regarding the GAS models and maximum likelihood estimation, see Blasques et al. (2014), Blasques et al. (2018), and Blasques et al. (2022).
#'
#' The use of the score for updating time-varying parameters is optimal in an information theoretic sense.
#' For an investigation of the optimality properties of GAS models, see Blasques et al. (2015) and Blasques et al. (2021).
#'
#' Generally, the GAS models perform quite well when compared to alternatives, including parameter-driven models.
#' For a comparison of the GAS models to alternative models, see Koopman et al. (2016) and Blazsek and Licht (2020).
#'
#' The GAS class includes many well-known econometric models, such as the generalized autoregressive conditional heteroskedasticity (GARCH) model of Bollerslev (1986), the autoregressive conditional duration (ACD) model of Engle and Russel (1998), and the Poisson count model of Davis et al. (2003).
#' More recently, a variety of novel score-driven models has been proposed, such as the Beta-t-(E)GARCH model of Harvey and Chakravarty (2008), the discrete price changes model of Koopman et al. (2018), the circular model of Harvey (2019), the bivariate Poisson model of Koopman and Lit (2019), and the ranking model of Holý and Zouhar (2022).
#' For an overview of various GAS models, see Harvey (2022).
#'
#' The extensive GAS literature is listed on \href{http://www.gasmodel.com}{www.gasmodel.com}.
#'
#' @return A \code{list} of S3 class \code{gas} with components:
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
#' \item{model$num_obs}{The actual number of observations used in the likelihood.}
#' \item{model$num_coef}{The actual number of estimated coefficients.}
#' \item{control$optim_function}{The optimization function.}
#' \item{control$optim_arguments}{The arguments which are passed to the optimization function.}
#' \item{control$hessian_function}{The Hessian function.}
#' \item{control$hessian_arguments}{The arguments which are passed to the Hessian function.}
#' \item{solution$status_start}{The status of the starting values computation.}
#' \item{solution$theta_start}{The computed starting values.}
#' \item{solution$status_optim}{The status of the optimization computation.}
#' \item{solution$theta_optim}{The computed optimal values.}
#' \item{solution$status_hessian}{The status of the Hessian computation.}
#' \item{solution$theta_hessian}{The computed Hessian.}
#' \item{fit$coef_est}{The estimated coefficients.}
#' \item{fit$coef_vcov}{The estimated variance-covariance matrix.}
#' \item{fit$coef_sd}{The estimated standard deviations.}
#' \item{fit$coef_zstat}{The statistics of the Z-test.}
#' \item{fit$coef_pval}{The p-values of the Z-test.}
#' \item{fit$par_unc}{The unconditional values of time-varying parameters.}
#' \item{fit$par_tv}{The individual values of time-varying parameter.}
#' \item{fit$score_tv}{The individual scores of time-varying parameters.}
#' \item{fit$mean_tv}{The expected values given by the model.}
#' \item{fit$var_tv}{The variances given by the model.}
#' \item{fit$resid_tv}{The residuals of the model.}
#' \item{fit$loglik_tv}{The log-likelihoods for the individual observations.}
#' \item{fit$loglik_sum}{The overall log-likelihood.}
#' \item{fit$aic}{The Akaike information criterion.}
#' \item{fit$bic}{The Bayesian information criterion.}
#'
#' @note
#' Supported generic functions for S3 class \code{gas} include \code{\link[base:summary]{summary()}}, \code{\link[base:plot]{plot()}}, \code{\link[stats:coef]{coef()}}, \code{\link[stats:vcov]{vcov()}}, \code{\link[stats:fitted]{fitted()}}, \code{\link[stats:residuals]{residuals()}}, \code{\link[stats:logLik]{logLik()}}, \code{\link[stats:AIC]{AIC()}}, \code{\link[stats:BIC]{BIC()}}, and \code{\link[stats:confint]{confint()}}.
#'
#' @references
#' Blasques, F., Gorgi, P., Koopman, S. J., and Wintenberger, O. (2018). Feasible Invertibility Conditions and Maximum Likelihood Estimation for Observation-Driven Models. \emph{Electronic Journal of Statistics}, \strong{12}(1), 1019–1052. \doi{10.1214/18-ejs1416}.
#'
#' Blasques, F., Koopman, S. J., and Lucas, A. (2014). Stationarity and Ergodicity of Univariate Generalized Autoregressive Score Processes. \emph{Electronic Journal of Statistics}, \strong{8}(1), 1088–1112. \doi{10.1214/14-ejs924}.
#'
#' Blasques, F., Koopman, S. J., and Lucas, A. (2015). Information-Theoretic Optimality of Observation-Driven Time Series Models for Continuous Responses. \emph{Biometrika}, \strong{102}(2), 325–343. \doi{10.1093/biomet/asu076}.
#'
#' Blasques, F., Lucas, A., and van Vlodrop, A. C. (2021). Finite Sample Optimality of Score-Driven Volatility Models: Some Monte Carlo Evidence. \emph{Econometrics and Statistics}, \strong{19}, 47–57. \doi{10.1016/j.ecosta.2020.03.010}.
#'
#' Blasques, F., van Brummelen, J., Koopman, S. J., and Lucas, A. (2022). Maximum Likelihood Estimation for Score-Driven Models. \emph{Journal of Econometrics}, \strong{227}(2), 325–346. \doi{10.1016/j.jeconom.2021.06.003}.
#'
#' Blazsek, S. and Licht, A. (2020). Dynamic Conditional Score Models: A Review of Their Applications. \emph{Applied Economics}, \strong{52}(11), 1181–1199. \doi{10.1080/00036846.2019.1659498}.
#'
#' Bollerslev, T. (1986). Generalized Autoregressive Conditional Heteroskedasticity. \emph{Journal of Econometrics}, \strong{31}(3), 307–327. \doi{10.1016/0304-4076(86)90063-1}.
#'
#' Creal, D., Koopman, S. J., and Lucas, A. (2013). Generalized Autoregressive Score Models with Applications. \emph{Journal of Applied Econometrics}, \strong{28}(5), 777–795. \doi{10.1002/jae.1279}.
#'
#' Davis, R. A., Dunsmuir, W. T. M., and Street, S. B. (2003). Observation-Driven Models for Poisson Counts. \emph{Biometrika}, \strong{90}(4), 777–790. \doi{10.1093/biomet/90.4.777}.
#'
#' Engle, R. F. and Russell, J. R. (1998). Autoregressive Conditional Duration: A New Model for Irregularly Spaced Transaction Data. \emph{Econometrica}, \strong{66}(5), 1127–1162. \doi{10.2307/2999632}.
#'
#' Harvey, A. C. (2013). \emph{Dynamic Models for Volatility and Heavy Tails: With Applications to Financial and Economic Time Series}. Cambridge University Press. \doi{10.1017/cbo9781139540933}.
#'
#' Harvey, A. C. (2022). Score-Driven Time Series Models. \emph{Annual Review of Statistics and Its Application}, \strong{9}(1), 321–342. \doi{10.1146/annurev-statistics-040120-021023}.
#'
#' Harvey, A. C. and Chakravarty, T. (2008). Beta-t-(E)GARCH. \emph{Cambridge Working Papers in Economics}, CWPE 0840. \doi{10.17863/cam.5286}.
#'
#' Harvey, A., Hurn, S., and Thiele, S. (2019). Modeling Directional (Circular) Time Series. \emph{Cambridge Working Papers in Economics}, CWPE 1971. \doi{10.17863/cam.43915}.
#'
#' Holý, V. and Zouhar, J. (2022). Modelling Time-Varying Rankings with Autoregressive and Score-Driven Dynamics. Journal of the Royal Statistical Society: Series C (Applied Statistics), \strong{71}(5). \doi{10.1111/rssc.12584}.
#'
#' Koopman, S. J. and Lit, R. (2019). Forecasting Football Match Results in National League Competitions Using Score-Driven Time Series Models. \emph{International Journal of Forecasting}, \strong{35}(2), 797–809. \doi{10.1016/j.ijforecast.2018.10.011}.
#'
#' Koopman, S. J., Lit, R., Lucas, A., and Opschoor, A. (2018). Dynamic Discrete Copula Models for High-Frequency Stock Price Changes. \emph{Journal of Applied Econometrics}, \strong{33}(7), 966–985. \doi{10.1002/jae.2645}.
#'
#' Koopman, S. J., Lucas, A., and Scharth, M. (2016). Predicting Time-Varying Parameters with Parameter-Driven and Observation-Driven Models. \emph{Review of Economics and Statistics}, \strong{98}(1), 97–110. \doi{10.1162/rest_a_00533}.
#'
#' @seealso
#' \code{\link[gasmodel:distr]{distr()}},
#' \code{\link[gasmodel:gas_bootstrap]{gas_bootstrap()}},
#' \code{\link[gasmodel:gas_filter]{gas_filter()}},
#' \code{\link[gasmodel:gas_forecast]{gas_forecast()}},
#' \code{\link[gasmodel:gas_simulate]{gas_simulate()}},
#' \code{\link[gasmodel:wrappers_optim]{wrappers_optim}},
#' \code{\link[gasmodel:wrappers_hessian]{wrappers_hessian}}
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
#' # Obtain the estimated coefficients
#' coef(est_negbin)
#'
#' # Obtain the estimated variance-covariance matrix
#' vcov(est_negbin)
#'
#' # Obtain the log-likelihood, AIC, and BIC
#' logLik(est_negbin)
#' AIC(est_negbin)
#' BIC(est_negbin)
#'
#' # Obtain the confidence intervals of coefficients
#' confint(est_negbin)
#'
#' # Plot the time-varying parameters
#' plot(est_negbin)}
#'
#' @export
gas <- function(y, x = NULL, distr, param = NULL, scaling = "unit", regress = "joint", p = 1L, q = 1L, par_static = NULL, par_link = NULL, par_init = NULL, lik_skip = 0L, coef_fix_value = NULL, coef_fix_other = NULL, coef_fix_special = NULL, coef_bound_lower = NULL, coef_bound_upper = NULL, coef_start = NULL, optim_function = wrapper_optim_nloptr, optim_arguments = list(opts = list(algorithm = 'NLOPT_LN_NELDERMEAD', xtol_rel = 0, maxeval = 1e6)), hessian_function = wrapper_hessian_stats, hessian_arguments = list(), print_progress = FALSE) {
  model <- list()
  model$distr <- check_my_distr(distr = distr)
  model$param <- check_my_param(param = param, distr = model$distr)
  model$scaling <- check_my_scaling(scaling = scaling)
  model$regress <- check_my_regress(regress = regress)
  info_distr <- info_distribution(distr = model$distr, param = model$param)
  data <- list()
  data$y <- check_my_y(y = y, dim = info_distr$dim, type = info_distr$type)
  model$t <- check_my_t(y = data$y)
  model$n <- check_my_n(y = data$y)
  info_par <- info_parameters(distr = model$distr, param = model$param, n = model$n)
  data$x <- check_my_x(x = x, t = model$t, par_num = info_par$par_num, group_num = info_par$group_num, par_in_group_num = info_par$par_in_group_num)
  model$m <- check_my_m(x = data$x)
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
  model$lik_skip <- check_my_lik_skip(lik_skip = lik_skip, t = model$t, p = model$p, q = model$q)
  info_coef <- info_coefficients(m = model$m, p = model$p, q = model$q, par_static = model$par_static, par_names = info_par$par_names, par_num = info_par$par_num, group_names = info_par$group_names, group_of_par_names = info_par$group_of_par_names)
  model$coef_fix_value <- check_my_coef_fix_value(coef_fix_value = coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_other <- check_my_coef_fix_other(coef_fix_other = coef_fix_other, coef_fix_value = model$coef_fix_value, coef_num = info_coef$coef_num)
  model$coef_fix_special <- check_my_coef_fix_special(coef_fix_special = coef_fix_special)
  model[c("coef_fix_value", "coef_fix_other")] <- fixed_coefficients(model = model, info_par = info_par, info_coef = info_coef)
  model$coef_bound_lower <- name_vector(check_my_coef_bound_lower(coef_bound_lower = coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  model$coef_bound_upper <- name_vector(check_my_coef_bound_upper(coef_bound_upper = coef_bound_upper, coef_bound_lower = model$coef_bound_lower, par_static = model$par_static, par_support = info_par$par_support, par_num = info_par$par_num, coef_in_par_num = info_coef$coef_in_par_num, coef_num = info_coef$coef_num), info_coef$coef_names)
  info_theta <- info_thetas(coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other, coef_names = info_coef$coef_names)
  fun <- list()
  fun$loglik <- setup_fun_loglik(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$mean <- setup_fun_mean(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$var <- setup_fun_var(distr = model$distr, param = model$param, par_trans = info_par$par_trans)
  fun$score <- setup_fun_score(distr = model$distr, param = model$param, scaling = model$scaling, orthog = info_distr$orthog, par_trans = info_par$par_trans, par_static = model$par_static)
  comp <- list()
  comp$coef_start <- check_my_coef_start(coef_start = coef_start, coef_bound_lower = model$coef_bound_lower, coef_bound_upper = model$coef_bound_upper, coef_num = info_coef$coef_num)
  comp$theta_start <- convert_coef_vector_to_theta_vector(comp$coef_start, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$theta_bound_lower <- convert_coef_vector_to_theta_vector(model$coef_bound_lower, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$theta_bound_upper <- convert_coef_vector_to_theta_vector(model$coef_bound_upper, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other)
  comp$compute_start <- any(is.na(comp$theta_start))
  control <- list()
  if (is.null(optim_function)) {
    comp$compute_optim <- FALSE
    control['optim_function'] <- list(NULL)
    control['optim_arguments'] <- list(NULL)
  } else {
    comp$compute_optim <- TRUE
    control$optim_function <- check_generic_function(arg = optim_function, arg_name = "optim_function")
    control$optim_arguments <- check_generic_list(arg = optim_arguments, arg_name = "optim_arguments")
  }
  if (is.null(hessian_function)) {
    comp$compute_hessian <- FALSE
    control['hessian_function'] <- list(NULL)
    control['hessian_arguments'] <- list(NULL)
  } else {
    comp$compute_hessian <- TRUE
    control$hessian_function <- check_generic_function(arg = hessian_function, arg_name = "hessian_function")
    control$hessian_arguments <- check_generic_list(arg = hessian_arguments, arg_name = "hessian_arguments")
  }
  comp$print_progress <- check_generic_logical_scalar(arg = print_progress, arg_name = "print_progress")
  comp$est_details <- list(data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, print_progress = comp$print_progress)
  solution <- list()
  if (comp$compute_start) {
    if (comp$print_progress) { message("Computing a starting solution...") }
    comp$result_start <- starting_theta(theta_start = solution$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, data = data, model = model, fun = fun, info_distr = info_distr, info_par = info_par, info_coef = info_coef, info_theta = info_theta, print_progress = comp$print_progress)
    solution$status_start <- comp$result_start$status_start
    solution$theta_start <- name_vector(comp$result_start$theta_start, info_theta$theta_names)
    if (solution$status_start != "success") {
      warning("Computation of a starting solution ended with status '", solution$status_start, "'.")
    }
  } else {
    solution$status_start <- "starting_values_supplied"
    solution$theta_start <- name_vector(comp$theta_start, info_theta$theta_names)
  }
  if (comp$compute_optim) {
    if (comp$print_progress) { message("Computing the optimal solution...") }
    comp$result_optim <- do.call(control$optim_function, args = c(list(obj_fun = likelihood_objective, theta_start = solution$theta_start, theta_bound_lower = comp$theta_bound_lower, theta_bound_upper = comp$theta_bound_upper, est_details = comp$est_details), control$optim_arguments))
    solution$status_optim <- comp$result_optim$status_optim
    solution$theta_optim <- name_vector(comp$result_optim$theta_optim, info_theta$theta_names)
    if (!(solution$status_optim %in% c("success", "objective_tolerance_reached", "desired_objective_reached", "variables_tolerance_reached"))) {
      warning("Computation of the optimal solution ended with status '", solution$status_optim, "'.")
    }
  } else {
    solution$status_optim <- "computation_skipped"
    solution$theta_optim <- solution$theta_start
  }
  if (comp$compute_hessian) {
    if (comp$print_progress) { message("Computing the Hessian matrix...") }
    comp$result_hessian <- do.call(control$hessian_function, args = c(list(obj_fun = likelihood_objective, theta_optim = solution$theta_optim, est_details = comp$est_details), control$hessian_arguments))
    solution$status_hessian <- comp$result_hessian$status_hessian
    solution$theta_hessian <- name_matrix(comp$result_hessian$theta_hessian, info_theta$theta_names, info_theta$theta_names)
    if (solution$status_hessian != "success") {
      warning("Computation of the Hessian matrix ended with status '", solution$status_hessian, "'.")
    }
  } else {
    solution$status_hessian <- "computation_skipped"
    solution$theta_hessian <- name_matrix(matrix(NA_real_, nrow = info_theta$theta_num, ncol = info_theta$theta_num), info_theta$theta_names, info_theta$theta_names)
  }
  fit <- list()
  fit$coef_est <- name_vector(convert_theta_vector_to_coef_vector(solution$theta_optim, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other), info_coef$coef_names)
  comp$eval_tv <- be_silent(likelihood_evaluate(coef = fit$coef_est, data = data, model = model, fun = fun, info_par = info_par, info_coef = info_coef))
  model$num_obs <- sum(!is.na(comp$eval_tv$lik))
  model$num_coef <- info_theta$theta_num
  comp$theta_vcov <- matrix_inv(solution$theta_hessian) / model$num_obs
  comp$struc <- convert_coef_vector_to_struc_list(coef_vec = fit$coef_est, m = model$m, p = model$p, q = model$q, par_names = info_par$par_names, par_of_coef_names = info_coef$par_of_coef_names)
  fit$coef_vcov <- name_matrix(convert_theta_matrix_to_coef_matrix(comp$theta_vcov, coef_fix_value = model$coef_fix_value, coef_fix_other = model$coef_fix_other), info_coef$coef_names, info_coef$coef_names)
  fit$coef_sd <- be_silent(sqrt(diag(fit$coef_vcov)))
  fit$coef_zstat <- fit$coef_est / fit$coef_sd
  fit$coef_pval <- 2 * stats::pnorm(-abs(fit$coef_zstat))
  if (model$regress == "joint") {
    fit$par_unc <- sapply(1:info_par$par_num, function(i) { (comp$struc[[i]]$omega + colMeans(data$x[[i]], na.rm = TRUE) %*% comp$struc[[i]]$beta) / (1 - sum(comp$struc[[i]]$phi)) })
  } else if (model$regress == "sep") {
    fit$par_unc <- sapply(1:info_par$par_num, function(i) { comp$struc[[i]]$omega + colMeans(data$x[[i]], na.rm = TRUE) %*% comp$struc[[i]]$beta })
  }
  info_data <- info_data(y = data$y, x = data$x)
  data$y <- name_matrix(data$y, info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
  data$x <- name_list_of_matrices(data$x, info_par$par_names, info_data$index_time_list, info_data$index_vars_list, drop = c(FALSE, TRUE), zero = c(FALSE, TRUE))
  fit$par_tv <- name_matrix(comp$eval_tv$par, info_data$index_time, info_par$par_names, drop = c(FALSE, TRUE))
  fit$score_tv <- name_matrix(comp$eval_tv$score, info_data$index_time, info_par$par_names, drop = c(FALSE, TRUE))
  fit$mean_tv <- name_matrix(fun$mean(comp$eval_tv$par), info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
  fit$var_tv <- name_matrix(convert_varcov_array_to_var_matrix(fun$var(comp$eval_tv$par)), info_data$index_time, info_data$index_series, drop = c(FALSE, TRUE))
  fit$resid_tv <- (data$y - fit$mean_tv) / sqrt(fit$var_tv)
  fit$loglik_tv <- name_vector(comp$eval_tv$lik, info_data$index_time)
  fit$loglik_sum <- sum(fit$loglik_tv, na.rm = TRUE)
  fit$aic <- 2 * model$num_coef - 2 * fit$loglik_sum
  fit$bic <- log(model$num_obs) * model$num_coef - 2 * fit$loglik_sum
  if (mean(fit$loglik_tv, na.rm = TRUE) <= -1e100) {
    warning("The likelihood function has zero value. The results are not reliable.")
  }
  report <- list(data = data, model = model, control = control, solution = solution, fit = fit)
  class(report) <- "gas"
  return(report)
}
# ------------------------------------------------------------------------------


# Print Estimate ---------------------------------------------------------------
#' @export
print.gas <- function(x, ...) {
  info_title <- info_title(distr = x$model$distr, param = x$model$param, scaling = x$model$scaling)
  cat("GAS Model:", info_title$title, "\n")
  cat("\n")
  cat("Coefficients:", "\n")
  coef_table <- cbind("Estimate" = x$fit$coef_est, "Std. Error" = x$fit$coef_sd, "Z-Test" = x$fit$coef_zstat, "Pr(>|Z|)" = x$fit$coef_pval)
  stats::printCoefmat(coef_table)
  cat("\n")
  cat("Log-Likelihood:", format(x$fit$loglik_sum))
  cat(", AIC:", format(x$fit$aic))
  cat(", BIC:", format(x$fit$bic))
  cat("\n")
  invisible(x)
}
# ------------------------------------------------------------------------------


# Summarize Estimate -----------------------------------------------------------
#' @export
summary.gas <- function(object, ...) {
  print(object)
  cat("\n")
  cat("Unconditional Parameters:", "\n")
  print(object$fit$par_unc)
  cat("\n")
  cat("Time-Varying Parameters:", "\n")
  print(object$fit$par_tv)
  invisible(object)
}
# ------------------------------------------------------------------------------


# Obtain Coefficients ----------------------------------------------------------
#' @export
coef.gas <- function(object, ...) {
  coef_est <- object$fit$coef_est
  return(coef_est)
}
# ------------------------------------------------------------------------------


# Obtain Variance-Covariance Matrix --------------------------------------------
#' @export
vcov.gas <- function(object, ...) {
  coef_vcov <- object$fit$coef_vcov
  return(coef_vcov)
}
# ------------------------------------------------------------------------------


# Obtain Fitted Values ---------------------------------------------------------
#' @export
fitted.gas <- function(object, ...) {
  mean_tv <- object$fit$mean_tv
  return(mean_tv)
}
# ------------------------------------------------------------------------------


# Obtain Residuals -------------------------------------------------------------
#' @export
residuals.gas <- function(object, ...) {
  resid <- object$fit$resid
  return(resid)
}
# ------------------------------------------------------------------------------


# Obtain Log-Likelihood --------------------------------------------------------
#' @export
logLik.gas <- function(object, ...) {
  loglik_sum <- object$fit$loglik_sum
  attr(loglik_sum, "nobs") <- object$model$num_obs
  attr(loglik_sum, "df") <- object$model$num_coef
  class(loglik_sum) <- "logLik"
  return(loglik_sum)
}
# ------------------------------------------------------------------------------


# Plot Time-Varying Parameters -------------------------------------------------
#' @importFrom ggplot2 .data
#' @export
plot.gas <- function(x, which = NULL, ...) {
  par_static <- x$model$par_static
  par_tv <- as.matrix(x$fit$par_tv)
  par_unc <- x$fit$par_unc
  par_names <- names(par_unc)
  par_num <- length(par_unc)
  be_silent(ts_index <- as.numeric(rownames(par_tv)))
  if (any(is.na(ts_index)) || any(diff(ts_index) <= 0)) {
    ts_index <- 1:nrow(par_tv)
  }
  gg_list <- list()
  for (i in which(!par_static)) {
    gg_data <- data.frame(index = ts_index, value = par_tv[, i])
    gg_fig <- ggplot2::ggplot(gg_data, mapping = ggplot2::aes(.data$index, .data$value)) +
      ggplot2::geom_hline(yintercept = par_unc[i], linetype = "dashed") +
      ggplot2::geom_line(color = "#800000") +
      ggplot2::geom_point(color = "#800000") +
      ggplot2::labs(title = paste("Time-Varying Parameter", par_names[i]), x = "Time Index", y = "Parameter Value")
    gg_list <- append(gg_list, list(gg_fig))
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


