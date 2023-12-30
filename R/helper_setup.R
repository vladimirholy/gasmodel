
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
setup_fun_score <- function(distr, param, scaling, orthog, par_trans, par_static) {
  if (distr == "exp" && param == "scale" && all(par_trans == c("logarithmic"))) {
    fun <- function(y, f) { y[1, 1] * exp(-f[1, 1]) - 1 }
  } else if (distr == "gamma" && param == "scale" && scaling == "unit" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(y[1, 1] * exp(-f[1, 1]) - f[1, 2], log(y[1, 1]) - f[1, 1] - digamma(f[1, 2])) }
  } else if (distr == "gamma" && param == "scale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] - f[1, 2] * exp(f[1, 1])) / (sqrt(f[1, 2]) * exp(f[1, 1])), (log(y[1, 1]) - f[1, 1] - digamma(f[1, 2])) / sqrt(trigamma(f[1, 2]))) }
  } else if (distr == "gamma" && param == "scale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(y[1, 1] / f[1, 2] * exp(-f[1, 1]) - 1, (log(y[1, 1]) - f[1, 1] - digamma(f[1, 2])) / trigamma(f[1, 2])) }
  } else if (distr == "gengamma" && param == "scale" && scaling == "unit" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(f[1, 3] * ((y[1, 1] * exp(-f[1, 1]))^f[1, 3] - f[1, 2]), f[1, 3] * (log(y[1, 1]) - f[1, 1]) - digamma(f[1, 2]), (f[1, 2] - (y[1, 1] * exp(-f[1, 1]))^f[1, 3]) * (log(y[1, 1]) - f[1, 1]) + 1 / f[1, 3]) }
  } else if (distr == "gengamma" && param == "scale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("logarithmic", "identity", "identity"))) {
    fun <- function(y, f) { c(((y[1, 1] * exp(-f[1, 1]))^f[1, 3] - f[1, 2]) / sqrt(f[1, 2]), (f[1, 3] * (log(y[1, 1]) - f[1, 1]) - digamma(f[1, 2])) / sqrt(trigamma(f[1, 2])), (((f[1, 2] - (y[1, 1] * exp(-f[1, 1]))^f[1, 3]) * (log(y[1, 1]) - f[1, 1])) * f[1, 3] + 1) / sqrt(f[1, 2] * digamma(f[1, 2])^2 + 2 * digamma(f[1, 2]) + f[1, 2] * trigamma(f[1, 2]) + 1)) }
  } else if (distr == "gengamma" && param == "scale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("logarithmic", "identity", "identity"))) {
    fun <- function(y, f) { c(((y[1, 1] * exp(-f[1, 1]))^f[1, 3] - f[1, 2]) / (f[1, 2] * f[1, 3]), (f[1, 3] * (log(y[1, 1]) - f[1, 1]) - digamma(f[1, 2])) / trigamma(f[1, 2]), (((f[1, 2] - (y[1, 1] * exp(-f[1, 1]))^f[1, 3]) * (log(y[1, 1]) - f[1, 1])) * f[1, 3]^2 + f[1, 3]) / (f[1, 2] * digamma(f[1, 2])^2 + 2 * digamma(f[1, 2]) + f[1, 2] * trigamma(f[1, 2]) + 1)) }
  } else if (distr == "weibull" && param == "scale" && scaling == "unit" && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(f[1, 2] * ((y[1, 1] * exp(-f[1, 1]))^f[1, 2] - 1), (1 - (y[1, 1] * exp(-f[1, 1]))^f[1, 2]) * (log(y[1, 1]) - f[1, 1]) + 1 / f[1, 2]) }
  } else if (distr == "weibull" && param == "scale" && ((scaling == "fisher_inv_sqrt" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv_sqrt") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c((y[1, 1] * exp(-f[1, 1]))^f[1, 2] - 1, (((1 - (y[1, 1] * exp(-f[1, 1]))^f[1, 2]) * (log(y[1, 1]) - 1)) * f[1, 2] + 1) / sqrt(digamma(1)^2 + 2 * digamma(1) + trigamma(1) + 1)) }
  } else if (distr == "weibull" && param == "scale" && ((scaling == "fisher_inv" && sum(!par_static) <= 1L) || scaling == "diag_fisher_inv") && all(par_trans == c("logarithmic", "identity"))) {
    fun <- function(y, f) { c(((y[1, 1] * exp(-f[1, 1]))^f[1, 2] - 1) / (f[1, 2]), (((1 - (y[1, 1] * exp(-f[1, 1]))^f[1, 2]) * (log(y[1, 1]) - f[1, 1])) * f[1, 2]^2 + f[1, 2]) / (digamma(1)^2 + 2 * digamma(1) + trigamma(1) + 1)) }
  } else if (scaling == "unit") {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob_inv <- reparam_link_jacob_inv_str(par_trans)
    str_score <- paste0(str_jacob_inv, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if ((scaling == "fisher_inv_sqrt" && (orthog == TRUE || sum(!par_static) <= 1L)) || (scaling == "full_fisher_inv_sqrt" && orthog == TRUE) || scaling == "diag_fisher_inv_sqrt") {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_fisher_inv_sqrt_orig <- paste0("be_silent(diag(1 / sqrt(diag(as.matrix(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ]))), nrow = ", length(par_trans), "))")
    str_score <- paste0(str_fisher_inv_sqrt_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "fisher_inv_sqrt" && orthog == FALSE && sum(!par_static) >= 2L) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob_inv <- reparam_link_jacob_inv_str(par_trans)
    str_fisher_inv_sqrt_tilde <- paste0("matrix_part_inv_sqrt(", str_jacob_inv, " %*% ", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ] %*% ", str_jacob_inv, ", idx = ", par_static, ")")
    str_score <- paste0(str_fisher_inv_sqrt_tilde, " %*% ", str_jacob_inv, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "full_fisher_inv_sqrt" && orthog == FALSE) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob_inv <- reparam_link_jacob_inv_str(par_trans)
    str_fisher_inv_sqrt_tilde <- paste0("matrix_inv_sqrt(", str_jacob_inv, " %*% ", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ] %*% ", str_jacob_inv, ")")
    str_score <- paste0(str_fisher_inv_sqrt_tilde, " %*% ", str_jacob_inv, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if ((scaling == "fisher_inv" && (orthog == TRUE || sum(!par_static) <= 1L)) || (scaling == "full_fisher_inv" && orthog == TRUE) || scaling == "diag_fisher_inv") {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob <- reparam_link_jacob_str(par_trans)
    str_fisher_inv_orig <- paste0("be_silent(diag(1 / diag(as.matrix(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ])), nrow = ", length(par_trans), "))")
    str_score <- paste0(str_jacob, " %*% ", str_fisher_inv_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "fisher_inv" && orthog == FALSE && sum(!par_static) >= 2L) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob <- reparam_link_jacob_str(par_trans)
    str_fisher_inv_orig <- paste0("matrix_part_inv(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ], idx = ", par_static, ")")
    str_score <- paste0(str_jacob, " %*% ", str_fisher_inv_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
    str_fun <- paste0("function(y, f) { ", str_score, "}")
    fun <- eval(parse(text = str_fun))
  } else if (scaling == "full_fisher_inv" && orthog == FALSE) {
    str_f_orig <- reparam_tilde_to_orig_str(par_trans)
    str_jacob <- reparam_link_jacob_str(par_trans)
    str_fisher_inv_orig <- paste0("matrix_inv(", paste("distr", distr, param, "fisher", sep = "_"), "(f = ", str_f_orig, ")[1, , ])")
    str_score <- paste0(str_jacob, " %*% ", str_fisher_inv_orig, " %*% ", paste("distr", distr, param, "score", sep = "_"), "(y = y, f = ", str_f_orig, ")[1, ]")
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


