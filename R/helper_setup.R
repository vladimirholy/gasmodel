
# INTERNAL DISTRIBUTION SETUP FUNCTIONS


# Create the Log-Likelihood Function -------------------------------------------
setup_fun_loglik <- function(distr, param, par_trans) {
  fun <- function(y, f, trans = NULL) {
    f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
    res_loglik <- do.call(paste("distr", distr, param, "loglik", sep = "_"), args = list(y = y, f = f_orig))
    res_loglik <- as.vector(res_loglik)
    return(res_loglik)
  }
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Mean Function -----------------------------------------------------
setup_fun_mean <- function(distr, param, par_trans) {
  fun <- function(f) {
    f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
    res_mean <- do.call(paste("distr", distr, param, "mean", sep = "_"), args = list(f = f_orig))
    return(res_mean)
  }
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Variance Function -------------------------------------------------
setup_fun_var <- function(distr, param, par_trans) {
  fun <- function(f) {
    f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
    res_var <- do.call(paste("distr", distr, param, "var", sep = "_"), args = list(f = f_orig))
    return(res_var)
  }
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Score Function ----------------------------------------------------
setup_fun_score <- function(distr, param, scaling, orthog, par_trans) {
  if (scaling == "unit") {
    fun <- function(y, f) {
      f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
      jacob_inv <- reparam_link_jacob_inv(f_orig = f_orig, par_trans = par_trans)
      res_score <- t(jacob_inv) %*% do.call(paste("distr", distr, param, "score", sep = "_"), args = list(y = y, f = f_orig))[1, ]
      return(res_score)
    }
  } else if (scaling == "fisher_inv" && orthog == TRUE) {
    fun <- function(y, f) {
      f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
      jacob <- reparam_link_jacob(f_orig = f_orig, par_trans = par_trans)
      fisher_inv_orig <- matrix_diag_inv(as.matrix(do.call(paste("distr", distr, param, "fisher", sep = "_"), args = list(f = f_orig))[1, , ]))
      res_score <- jacob %*% fisher_inv_orig %*% do.call(paste("distr", distr, param, "score", sep = "_"), args = list(y = y, f = f_orig))[1, ]
      return(res_score)
    }
  } else if (scaling == "fisher_inv" && orthog == FALSE) {
    fun <- function(y, f) {
      f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
      jacob <- reparam_link_jacob(f_orig = f_orig, par_trans = par_trans)
      fisher_inv_orig <- matrix_inv(as.matrix(do.call(paste("distr", distr, param, "fisher", sep = "_"), args = list(f = f_orig))[1, , ]))
      res_score <- jacob %*% fisher_inv_orig %*% do.call(paste("distr", distr, param, "score", sep = "_"), args = list(y = y, f = f_orig))[1, ]
      return(res_score)
    }
  } else if (scaling == "fisher_inv_sqrt" && orthog == TRUE) {
    fun <- function(y, f) {
      f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
      jacob_inv <- reparam_link_jacob_inv(f_orig = f_orig, par_trans = par_trans)
      fisher_inv_sqrt_tilde <- matrix_diag_inv_sqrt(t(jacob_inv) %*% as.matrix(do.call(paste("distr", distr, param, "fisher", sep = "_"), args = list(f = f_orig))[1, , ]) %*% jacob_inv)
      res_score <- fisher_inv_sqrt_tilde %*% t(jacob_inv) %*% do.call(paste("distr", distr, param, "score", sep = "_"), args = list(y = y, f = f_orig))[1, ]
      return(res_score)
    }
  } else if (scaling == "fisher_inv_sqrt" && orthog == FALSE) {
    fun <- function(y, f) {
      f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
      jacob_inv <- reparam_link_jacob_inv(f_orig = f_orig, par_trans = par_trans)
      fisher_inv_sqrt_tilde <- matrix_inv_sqrt(t(jacob_inv) %*% as.matrix(do.call(paste("distr", distr, param, "fisher", sep = "_"), args = list(f = f_orig))[1, , ]) %*% jacob_inv)
      res_score <- fisher_inv_sqrt_tilde %*% t(jacob_inv) %*% do.call(paste("distr", distr, param, "score", sep = "_"), args = list(y = y, f = f_orig))[1, ]
      return(res_score)
    }
  }
  return(fun)
}
# ------------------------------------------------------------------------------


# Create the Random Generation Function ----------------------------------------
setup_fun_random <- function(distr, param, par_trans) {
  fun <- function(t, f) {
    f_orig <- reparam_tilde_to_orig(f_tilde = f, par_trans = par_trans)
    res_random <- do.call(paste("distr", distr, param, "random", sep = "_"), args = list(t = t, f = f_orig))
    return(res_random)
  }
  return(fun)
}
# ------------------------------------------------------------------------------


