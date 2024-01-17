
# INTERNAL DISTRIBUTION SETUP FUNCTIONS


# Create the Log-Likelihood Function -------------------------------------------
setup_fun_loglik <- function(distr, param, par_trans) {
  str_f_orig <- reparam_tilde_to_orig_str(par_trans)
  str_loglik <- paste0("as.vector(", paste("distr", distr, param, "loglik", sep = "_"), "(y = y, f = ", str_f_orig, "))")
  str_fun <- paste0("function(y, f) { ", str_loglik, "}")
  fun <- eval(parse(text = str_fun))
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Mean Function -----------------------------------------------------
setup_fun_mean <- function(distr, param, par_trans) {
  str_f_orig <- reparam_tilde_to_orig_str(par_trans)
  str_mean <- paste0(paste("distr", distr, param, "mean", sep = "_"), "(f = ", str_f_orig, ")")
  str_fun <- paste0("function(f) { ", str_mean, "}")
  fun <- eval(parse(text = str_fun))
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Variance Function -------------------------------------------------
setup_fun_var <- function(distr, param, par_trans) {
  str_f_orig <- reparam_tilde_to_orig_str(par_trans)
  str_var <- paste0(paste("distr", distr, param, "var", sep = "_"), "(f = ", str_f_orig, ")")
  str_fun <- paste0("function(f) { ", str_var, "}")
  fun <- eval(parse(text = str_fun))
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Score Function ----------------------------------------------------
setup_fun_score <- function(distr, param, scaling, orthog, par_trans, par_static, fast = TRUE) {
  if (fast && distr == "alaplace" && param == "meanscale" && scaling == "unit" && all(par_trans == c("identity", "identity", "identity"))) {
    fun <- function(y, f) { c(sign(y[1, 1] - f[1, 1]) * f[1, 3]^sign(y[1, 1] - f[1, 1]) / f[1, 2], abs(y[1, 1] - f[1, 1]) * f[1, 3]^sign(y[1, 1] - f[1, 1]) / f[1, 2]^2 - 1 / f[1, 2], -(y[1, 1] - f[1, 1]) * f[1, 3]^(sign(y[1, 1] - f[1, 1]) - 1) / f[1, 2] + (1 - f[1, 3]^2) / (f[1, 3] * (1 + f[1, 3]^2))) }
  } else if (fast && distr == "alaplace" && param == "meanscale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("identity", "identity", "identity"))) {
    fun <- function(y, f) { c(f[1, 2] * (sign(y[1, 1] - f[1, 1]) * f[1, 3]^sign(y[1, 1] - f[1, 1]) / f[1, 2]), f[1, 2] * (abs(y[1, 1] - f[1, 1]) * f[1, 3]^sign(y[1, 1] - f[1, 1]) / f[1, 2]^2 - 1 / f[1, 2]), (-(y[1, 1] - f[1, 1]) * f[1, 3]^(sign(y[1, 1] - f[1, 1]) - 1) / f[1, 2] + (1 - f[1, 3]^2) / (f[1, 3] * (1 + f[1, 3]^2))) / sqrt((1 / f[1, 3]^2 + 4 / (1 + f[1, 3]^2)^2))) }
  } else if (fast && distr == "alaplace" && param == "meanscale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("identity", "identity", "identity"))) {
    fun <- function(y, f) { c(f[1, 2]^2 * (sign(y[1, 1] - f[1, 1]) * f[1, 3]^sign(y[1, 1] - f[1, 1]) / f[1, 2]), f[1, 2]^2 * (abs(y[1, 1] - f[1, 1]) * f[1, 3]^sign(y[1, 1] - f[1, 1]) / f[1, 2]^2 - 1 / f[1, 2]), (-(y[1, 1] - f[1, 1]) * f[1, 3]^(sign(y[1, 1] - f[1, 1]) - 1) / f[1, 2] + (1 - f[1, 3]^2) / (f[1, 3] * (1 + f[1, 3]^2))) / (1 / f[1, 3]^2 + 4 / (1 + f[1, 3]^2)^2)) }
  } else if (fast && distr == "bernoulli" && param == "prob" && scaling == "unit" && all(par_trans == c("logit"))) {
    fun <- function(y, f) { (y[1, 1] == 0) * (-1) / (1 + exp(-f[1, 1])) + (y[1, 1] == 1) / (1 + exp(f[1, 1])) }
  } else if (fast && distr == "bernoulli" && param == "prob" && scaling == "fisher_inv_sqrt" && all(par_trans == c("logit"))) {
    fun <- function(y, f) { (y[1, 1] == 0) * (-exp(f[1, 1] / 2)) + (y[1, 1] == 1) * (exp(-f[1, 1] / 2))  }
  } else if (fast && distr == "bernoulli" && param == "prob" && scaling == "fisher_inv" && all(par_trans == c("logit"))) {
    fun <- function(y, f) { (y[1, 1] == 0) * (-exp(f[1, 1]) - 1) + (y[1, 1] == 1) * (exp(-f[1, 1]) + 1) }
  } else if (fast && distr == "exp" && param == "scale" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { y[1, 1] * exp(-f[1, 1]) - 1 }
  } else if (fast && distr == "gamma" && param == "scale" && scaling == "unit" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(y[1, 1] * exp(-f[1, 1]) - f[1, 2], log(y[1, 1]) - f[1, 1] - digamma(f[1, 2])) }
  } else if (fast && distr == "gamma" && param == "scale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 2] * exp(f[1, 1])) / (sqrt(f[1, 2]) * exp(f[1, 1])), (log(y[1, 1]) - f[1, 1] - digamma(f[1, 2])) / sqrt(trigamma(f[1, 2]))) }
  } else if (fast && distr == "gamma" && param == "scale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(y[1, 1] / f[1, 2] * exp(-f[1, 1]) - 1, (log(y[1, 1]) - f[1, 1] - digamma(f[1, 2])) / trigamma(f[1, 2])) }
  } else if (fast && distr == "gengamma" && param == "scale" && scaling == "unit" && all(par_trans == c("logarithmic", "identity", "identity"))) {
    fun <- function(y, f) { c(f[1, 3] * ((y[1, 1] * exp(-f[1, 1]))^f[1, 3] - f[1, 2]), f[1, 3] * (log(y[1, 1]) - f[1, 1]) - digamma(f[1, 2]), (f[1, 2] - (y[1, 1] * exp(-f[1, 1]))^f[1, 3]) * (log(y[1, 1]) - f[1, 1]) + 1 / f[1, 3]) }
  } else if (fast && distr == "gengamma" && param == "scale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("logarithmic", "identity", "identity"))) {
    fun <- function(y, f) { c(((y[1, 1] * exp(-f[1, 1]))^f[1, 3] - f[1, 2]) / sqrt(f[1, 2]), (f[1, 3] * (log(y[1, 1]) - f[1, 1]) - digamma(f[1, 2])) / sqrt(trigamma(f[1, 2])), (((f[1, 2] - (y[1, 1] * exp(-f[1, 1]))^f[1, 3]) * (log(y[1, 1]) - f[1, 1])) * f[1, 3] + 1) / sqrt(f[1, 2] * digamma(f[1, 2])^2 + 2 * digamma(f[1, 2]) + f[1, 2] * trigamma(f[1, 2]) + 1)) }
  } else if (fast && distr == "gengamma" && param == "scale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("logarithmic", "identity", "identity"))) {
    fun <- function(y, f) { c(((y[1, 1] * exp(-f[1, 1]))^f[1, 3] - f[1, 2]) / (f[1, 2] * f[1, 3]), (f[1, 3] * (log(y[1, 1]) - f[1, 1]) - digamma(f[1, 2])) / trigamma(f[1, 2]), (((f[1, 2] - (y[1, 1] * exp(-f[1, 1]))^f[1, 3]) * (log(y[1, 1]) - f[1, 1])) * f[1, 3]^2 + f[1, 3]) / (f[1, 2] * digamma(f[1, 2])^2 + 2 * digamma(f[1, 2]) + f[1, 2] * trigamma(f[1, 2]) + 1)) }
  } else if (fast && distr == "geom" && param == "mean" && scaling == "unit" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { (y[1, 1] - exp(f[1, 1])) / (exp(f[1, 1]) + 1) }
  } else if (fast && distr == "geom" && param == "mean" && scaling == "fisher_inv_sqrt" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { (y[1, 1] - exp(f[1, 1])) / sqrt(exp(f[1, 1])^2 + exp(f[1, 1])) }
  } else if (fast && distr == "geom" && param == "mean" && scaling == "fisher_inv" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { y[1, 1] / exp(f[1, 1]) - 1 }
  } else if (fast && distr == "laplace" && param == "meanscale" && scaling == "unit" && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(sign(y[1, 1] - f[1, 1]) / f[1, 2], abs(y[1, 1] - f[1, 1]) / f[1, 2]^2 - 1 / f[1, 2]) }
  } else if (fast && distr == "laplace" && param == "meanscale" && (scaling == "fisher_inv_sqrt" || scaling == "diag_fisher_inv_sqrt" || scaling == "full_fisher_inv_sqrt") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(sign(y[1, 1] - f[1, 1]), abs(y[1, 1] - f[1, 1]) / f[1, 2] - 1) }
  } else if (fast && distr == "laplace" && param == "meanscale" && (scaling == "fisher_inv" || scaling == "diag_fisher_inv" || scaling == "full_fisher_inv") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(sign(y[1, 1] - f[1, 1]) * f[1, 2], abs(y[1, 1] - f[1, 1]) - f[1, 2]) }
  } else if (fast && distr == "logistic" && param == "meanscale" && scaling == "unit" && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(tanh((y[1, 1] - f[1, 1]) / (2 * f[1, 2])) / f[1, 2], (y[1, 1] - f[1, 1]) / f[1, 2]^2 * tanh((y[1, 1] - f[1, 1]) / (2 * f[1, 2])) - 1 / f[1, 2]) }
  } else if (fast && distr == "logistic" && param == "meanscale" && (scaling == "fisher_inv_sqrt" || scaling == "diag_fisher_inv_sqrt" || scaling == "full_fisher_inv_sqrt") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(tanh((y[1, 1] - f[1, 1]) / (2 * f[1, 2])) * sqrt(3), ((y[1, 1] - f[1, 1]) / f[1, 2] * tanh((y[1, 1] - f[1, 1]) / (2 * f[1, 2])) * sqrt(3) - sqrt(3)) / sqrt(pi^2 / 3 + 1)) }
  } else if (fast && distr == "logistic" && param == "meanscale" && (scaling == "fisher_inv" || scaling == "diag_fisher_inv" || scaling == "full_fisher_inv") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(tanh((y[1, 1] - f[1, 1]) / (2 * f[1, 2])) * 3 * f[1, 2], ((y[1, 1] - f[1, 1]) * tanh((y[1, 1] - f[1, 1]) / (2 * f[1, 2])) * 3 - 3 * f[1, 2]) / (pi^2 / 3 + 1)) }
  } else if (fast && distr == "negbin" && param == "nb2" && scaling == "unit" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - exp(f[1, 1])) / (1 + f[1, 2] * exp(f[1, 1])), (y[1, 1] - exp(f[1, 1])) / f[1, 2] / (1 + f[1, 2] * exp(f[1, 1])) + (log(1 + f[1, 2] * exp(f[1, 1])) + digamma(1 / f[1, 2]) - digamma(y + 1 / f[1, 2])) / f[1, 2]^2) }
  } else if (fast && distr == "negbin" && param == "nb2" && scaling == "fisher_inv_sqrt" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - exp(f[1, 1])) / exp(f[1, 1]) / (1 + f[1, 2] * exp(f[1, 1])) * sqrt(exp(f[1, 1]) * (1 + f[1, 2] * exp(f[1, 1]))), ((y[1, 1] - exp(f[1, 1])) / f[1, 2] / (1 + f[1, 2] * exp(f[1, 1])) + (log(1 + f[1, 2] * exp(f[1, 1])) + digamma(1 / f[1, 2]) - digamma(y + 1 / f[1, 2])) / f[1, 2]^2) / sqrt((log(1 + f[1, 2] * exp(f[1, 1])) + digamma(1 / f[1, 2]) - digamma(exp(f[1, 1]) + 1 / f[1, 2]))^2)) }
  } else if (fast && distr == "negbin" && param == "nb2" && scaling == "fisher_inv" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - exp(f[1, 1])) / exp(2 * f[1, 1]) / (1 + f[1, 2] * exp(f[1, 1])) * (exp(f[1, 1]) * (1 + f[1, 2] * exp(f[1, 1]))), ((y[1, 1] - exp(f[1, 1])) / f[1, 2] / (1 + f[1, 2] * exp(f[1, 1])) + (log(1 + f[1, 2] * exp(f[1, 1])) + digamma(1 / f[1, 2]) - digamma(y + 1 / f[1, 2])) / f[1, 2]^2) / ((log(1 + f[1, 2] * exp(f[1, 1])) + digamma(1 / f[1, 2]) - digamma(exp(f[1, 1]) + 1 / f[1, 2]))^2)) }
  } else if (fast && distr == "norm" && param == "meanvar" && scaling == "unit" && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 1]) / f[1, 2], ((y[1, 1] - f[1, 1])^2 - f[1, 2]) / (2 * f[1, 2]^2)) }
  } else if (fast && distr == "norm" && param == "meanvar" && (scaling == "fisher_inv_sqrt" || scaling == "diag_fisher_inv_sqrt" || scaling == "full_fisher_inv_sqrt") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 1]) / sqrt(f[1, 2]), ((y[1, 1] - f[1, 1])^2 - f[1, 2]) / (sqrt(2) * f[1, 2])) }
  } else if (fast && distr == "norm" && param == "meanvar" && (scaling == "fisher_inv" || scaling == "diag_fisher_inv" || scaling == "full_fisher_inv") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 1]), ((y[1, 1] - f[1, 1])^2 - f[1, 2])) }
  } else if (fast && distr == "norm" && param == "meanvar" && scaling == "unit" && all(par_trans == c("identity", "logarithmic"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 1]) / exp(f[1, 2]), ((y[1, 1] - f[1, 1])^2 - exp(f[1, 2])) / (2 * exp(f[1, 2]))) }
  } else if (fast && distr == "norm" && param == "meanvar" && (scaling == "fisher_inv_sqrt" || scaling == "diag_fisher_inv_sqrt" || scaling == "full_fisher_inv_sqrt") && all(par_trans == c("identity", "logarithmic"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 1]) / exp(f[1, 2] / 2), ((y[1, 1] - f[1, 1])^2 - exp(f[1, 2])) / (sqrt(2) * exp(f[1, 2]))) }
  } else if (fast && distr == "norm" && param == "meanvar" && (scaling == "fisher_inv" || scaling == "diag_fisher_inv" || scaling == "full_fisher_inv") && all(par_trans == c("identity", "logarithmic"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 1]), ((y[1, 1] - f[1, 1])^2 - exp(f[1, 2])) / exp(f[1, 2])) }
  } else if (fast && distr == "pois" && param == "mean" && scaling == "unit" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { y[1, 1] - exp(f[1, 1]) }
  } else if (fast && distr == "pois" && param == "mean" && scaling == "fisher_inv_sqrt" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { (y[1, 1] - exp(f[1, 1])) * exp(-f[1, 1] / 2) }
  } else if (fast && distr == "pois" && param == "mean" && scaling == "fisher_inv" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { y[1, 1] / exp(f[1, 1]) - 1 }
  } else if (fast && distr == "t" && param == "meanvar" && scaling == "unit" && all(par_trans == c("identity", "identity", "identity"))) {
    fun <- function(y, f) { c((f[1, 3] + 1) * (y[1, 1] - f[1, 1]) / (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * f[1, 2] + y[1, 1]^2), f[1, 3] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] - f[1, 2] + y[1, 1]^2) / (2 * f[1, 2] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * f[1, 2] + y[1, 1]^2)), ((y[1, 1] - f[1, 1])^2 - f[1, 2]) / ((y[1, 1] - f[1, 1])^2 + f[1, 2] * f[1, 3]) / 2 - log((y[1, 1] - f[1, 1])^2 / (f[1, 3] * f[1, 2]) + 1) / 2 - digamma(f[1, 3] / 2) / 2 + digamma((f[1, 3] + 1) / 2) / 2) }
  } else if (fast && distr == "t" && param == "meanvar" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("identity", "identity", "identity"))) {
    fun <- function(y, f) { c(((f[1, 3] + 1) * (y[1, 1] - f[1, 1]) / (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * f[1, 2] + y[1, 1]^2)) / sqrt((f[1, 3] + 1) / ((f[1, 3] + 3) * f[1, 2])), (f[1, 3] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] - f[1, 2] + y[1, 1]^2) / (2 * f[1, 2] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * f[1, 2] + y[1, 1]^2))) / sqrt(f[1, 3] / (2 * (f[1, 3] + 3) * f[1, 2]^2)), (((y[1, 1] - f[1, 1])^2 - f[1, 2]) / ((y[1, 1] - f[1, 1])^2 + f[1, 2] * f[1, 3]) / 2 - log((y[1, 1] - f[1, 1])^2 / (f[1, 3] * f[1, 2]) + 1) / 2 - digamma(f[1, 3] / 2) / 2 + digamma((f[1, 3] + 1) / 2) / 2) / sqrt(1 / 4 * trigamma(f[1, 3] / 2) - 1 / 4 * trigamma((f[1, 3] + 1) / 2) - (f[1, 3] + 5) / (2 * f[1, 3] * (f[1, 3] + 1) * (f[1, 3] + 3)))) }
  } else if (fast && distr == "t" && param == "meanvar" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("identity", "identity", "identity"))) {
    fun <- function(y, f) { c(((f[1, 3] + 1) * (y[1, 1] - f[1, 1]) / (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * f[1, 2] + y[1, 1]^2)) / ((f[1, 3] + 1) / ((f[1, 3] + 3) * f[1, 2])), (f[1, 3] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] - f[1, 2] + y[1, 1]^2) / (2 * f[1, 2] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * f[1, 2] + y[1, 1]^2))) / (f[1, 3] / (2 * (f[1, 3] + 3) * f[1, 2]^2)), (((y[1, 1] - f[1, 1])^2 - f[1, 2]) / ((y[1, 1] - f[1, 1])^2 + f[1, 2] * f[1, 3]) / 2 - log((y[1, 1] - f[1, 1])^2 / (f[1, 3] * f[1, 2]) + 1) / 2 - digamma(f[1, 3] / 2) / 2 + digamma((f[1, 3] + 1) / 2) / 2) / (1 / 4 * trigamma(f[1, 3] / 2) - 1 / 4 * trigamma((f[1, 3] + 1) / 2) - (f[1, 3] + 5) / (2 * f[1, 3] * (f[1, 3] + 1) * (f[1, 3] + 3)))) }
  } else if (fast && distr == "t" && param == "meanvar" && scaling == "unit" && all(par_trans == c("identity", "logarithmic", "identity"))) {
    fun <- function(y, f) { c((f[1, 3] + 1) * (y[1, 1] - f[1, 1]) / (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * exp(f[1, 2]) + y[1, 1]^2), f[1, 3] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] - exp(f[1, 2]) + y[1, 1]^2) / (2 * exp(f[1, 2]) * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * exp(f[1, 2]) + y[1, 1]^2)) * exp(f[1, 2]), ((y[1, 1] - f[1, 1])^2 - exp(f[1, 2])) / ((y[1, 1] - f[1, 1])^2 + exp(f[1, 2]) * f[1, 3]) / 2 - log((y[1, 1] - f[1, 1])^2 / (f[1, 3] * exp(f[1, 2])) + 1) / 2 - digamma(f[1, 3] / 2) / 2 + digamma((f[1, 3] + 1) / 2) / 2) }
  } else if (fast && distr == "t" && param == "meanvar" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("identity", "logarithmic", "identity"))) {
    fun <- function(y, f) { c(((f[1, 3] + 1) * (y[1, 1] - f[1, 1]) / (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * exp(f[1, 2]) + y[1, 1]^2)) / sqrt((f[1, 3] + 1) / ((f[1, 3] + 3) * exp(f[1, 2]))), (f[1, 3] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] - exp(f[1, 2]) + y[1, 1]^2) / (2 * exp(f[1, 2]) * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * exp(f[1, 2]) + y[1, 1]^2))) / sqrt(f[1, 3] / (2 * (f[1, 3] + 3) * exp(2 * f[1, 2]))), (((y[1, 1] - f[1, 1])^2 - exp(f[1, 2])) / ((y[1, 1] - f[1, 1])^2 + exp(f[1, 2]) * f[1, 3]) / 2 - log((y[1, 1] - f[1, 1])^2 / (f[1, 3] * exp(f[1, 2])) + 1) / 2 - digamma(f[1, 3] / 2) / 2 + digamma((f[1, 3] + 1) / 2) / 2) / sqrt(1 / 4 * trigamma(f[1, 3] / 2) - 1 / 4 * trigamma((f[1, 3] + 1) / 2) - (f[1, 3] + 5) / (2 * f[1, 3] * (f[1, 3] + 1) * (f[1, 3] + 3)))) }
  } else if (fast && distr == "t" && param == "meanvar" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("identity", "logarithmic", "identity"))) {
    fun <- function(y, f) { c(((f[1, 3] + 1) * (y[1, 1] - f[1, 1]) / (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * exp(f[1, 2]) + y[1, 1]^2)) / ((f[1, 3] + 1) / ((f[1, 3] + 3) * exp(f[1, 2]))), (f[1, 3] * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] - exp(f[1, 2]) + y[1, 1]^2) / (2 * exp(f[1, 2]) * (f[1, 1]^2 - 2 * f[1, 1] * y[1, 1] + f[1, 3] * exp(f[1, 2]) + y[1, 1]^2))) / (f[1, 3] / (2 * (f[1, 3] + 3) * exp(2 * f[1, 2]))) / exp(f[1, 2]), (((y[1, 1] - f[1, 1])^2 - exp(f[1, 2])) / ((y[1, 1] - f[1, 1])^2 + exp(f[1, 2]) * f[1, 3]) / 2 - log((y[1, 1] - f[1, 1])^2 / (f[1, 3] * exp(f[1, 2])) + 1) / 2 - digamma(f[1, 3] / 2) / 2 + digamma((f[1, 3] + 1) / 2) / 2) / (1 / 4 * trigamma(f[1, 3] / 2) - 1 / 4 * trigamma((f[1, 3] + 1) / 2) - (f[1, 3] + 5) / (2 * f[1, 3] * (f[1, 3] + 1) * (f[1, 3] + 3)))) }
  } else if (fast && distr == "vonmises" && param == "meanconc" && scaling == "unit" && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c(f[1, 2] * sin(y[1, 1] - f[1, 1]), cos(y[1, 1] - f[1, 1]) - besselI(f[1, 2], nu = 1) / besselI(f[1, 2], nu = 0)) }
  } else if (fast && distr == "vonmises" && param == "meanconc" && (scaling == "fisher_inv_sqrt" || scaling == "diag_fisher_inv_sqrt" || scaling == "full_fisher_inv_sqrt") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c((f[1, 2] * sin(y[1, 1] - f[1, 1])) / sqrt(f[1, 2] * besselI(f[1, 2], nu = 1) / besselI(f[1, 2], nu = 0)), (cos(y[1, 1] - f[1, 1]) - besselI(f[1, 2], nu = 1) / besselI(f[1, 2], nu = 0)) / sqrt(1 / 2 - besselI(f[1, 2], nu = 1)^2 / besselI(f[1, 2], nu = 0)^2 + besselI(f[1, 2], nu = 2) / (2 * besselI(f[1, 2], nu = 0)))) }
  } else if (fast && distr == "vonmises" && param == "meanconc" && (scaling == "fisher_inv" || scaling == "diag_fisher_inv" || scaling == "full_fisher_inv") && all(par_trans == c("identity", "identity"))) {
    fun <- function(y, f) { c((f[1, 2] * sin(y[1, 1] - f[1, 1])) / (f[1, 2] * besselI(f[1, 2], nu = 1) / besselI(f[1, 2], nu = 0)), (cos(y[1, 1] - f[1, 1]) - besselI(f[1, 2], nu = 1) / besselI(f[1, 2], nu = 0)) / (1 / 2 - besselI(f[1, 2], nu = 1)^2 / besselI(f[1, 2], nu = 0)^2 + besselI(f[1, 2], nu = 2) / (2 * besselI(f[1, 2], nu = 0)))) }
  } else if (fast && distr == "weibull" && param == "scale" && scaling == "unit" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(f[1, 2] * ((y[1, 1] * exp(-f[1, 1]))^f[1, 2] - 1), (1 - (y[1, 1] * exp(-f[1, 1]))^f[1, 2]) * (log(y[1, 1]) - f[1, 1]) + 1 / f[1, 2]) }
  } else if (fast && distr == "weibull" && param == "scale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] * exp(-f[1, 1]))^f[1, 2] - 1, (((1 - (y[1, 1] * exp(-f[1, 1]))^f[1, 2]) * (log(y[1, 1]) - f[1, 1])) * f[1, 2] + 1) / sqrt(digamma(1)^2 + 2 * digamma(1) + trigamma(1) + 1)) }
  } else if (fast && distr == "weibull" && param == "scale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(((y[1, 1] * exp(-f[1, 1]))^f[1, 2] - 1) / (f[1, 2]), (((1 - (y[1, 1] * exp(-f[1, 1]))^f[1, 2]) * (log(y[1, 1]) - f[1, 1])) * f[1, 2]^2 + f[1, 2]) / (digamma(1)^2 + 2 * digamma(1) + trigamma(1) + 1)) }
  } else if (scaling == "unit") {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob_inv <- reparam_link_jacob_inv_str(par_trans)
    str_score <- paste0("as.vector(", str_jacob_inv, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if ((scaling == "fisher_inv_sqrt" && (orthog == TRUE || sum(!par_static) <= 1L)) || (scaling == "full_fisher_inv_sqrt" && orthog == TRUE) || scaling == "diag_fisher_inv_sqrt") {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_fisher_inv_sqrt_orig <- paste0("be_silent(diag(1 / sqrt(diag(as.matrix(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ]))), nrow = ", length(par_trans), "))")
    str_score <- paste0("as.vector(", str_fisher_inv_sqrt_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "fisher_inv_sqrt" && orthog == FALSE && sum(!par_static) >= 2L) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob_inv <- reparam_link_jacob_inv_str(par_trans)
    str_fisher_inv_sqrt_tilde <- paste0("matrix_part_inv_sqrt(", str_jacob_inv, " %*% ", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ] %*% ", str_jacob_inv, ", idx = ", par_static, ")")
    str_score <- paste0("as.vector(", str_fisher_inv_sqrt_tilde, " %*% ", str_jacob_inv, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "full_fisher_inv_sqrt" && orthog == FALSE) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob_inv <- reparam_link_jacob_inv_str(par_trans)
    str_fisher_inv_sqrt_tilde <- paste0("matrix_inv_sqrt(", str_jacob_inv, " %*% ", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ] %*% ", str_jacob_inv, ")")
    str_score <- paste0("as.vector(", str_fisher_inv_sqrt_tilde, " %*% ", str_jacob_inv, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if ((scaling == "fisher_inv" && (orthog == TRUE || sum(!par_static) <= 1L)) || (scaling == "full_fisher_inv" && orthog == TRUE) || scaling == "diag_fisher_inv") {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob <- reparam_link_jacob_str(par_trans)
    str_fisher_inv_orig <- paste0("be_silent(diag(1 / diag(as.matrix(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ])), nrow = ", length(par_trans), "))")
    str_score <- paste0("as.vector(", str_jacob, " %*% ", str_fisher_inv_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "fisher_inv" && orthog == FALSE && sum(!par_static) >= 2L) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob <- reparam_link_jacob_str(par_trans)
    str_fisher_inv_orig <- paste0("matrix_part_inv(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ], idx = ", par_static, ")")
    str_score <- paste0("as.vector(", str_jacob, " %*% ", str_fisher_inv_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "full_fisher_inv" && orthog == FALSE) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob <- reparam_link_jacob_str(par_trans)
    str_fisher_inv_orig <- paste0("matrix_inv(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ])")
    str_score <- paste0("as.vector(", str_jacob, " %*% ", str_fisher_inv_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ])")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  }
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Random Generation Function ----------------------------------------
setup_fun_random <- function(distr, param, par_trans) {
  str_f_orig <- reparam_tilde_to_orig_str(par_trans)
  str_random <- paste0(paste("distr", distr, param, "random", sep = "_"), "(t = t, f = ", str_f_orig, ")")
  str_fun <- paste0("function(t, f) { ", str_random, "}")
  fun <- eval(parse(text = str_fun))
  return(fun)
}
# ------------------------------------------------------------------------------


