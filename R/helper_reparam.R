
# REPARAMETRIZATION FUNCTIONS


# Convert Tilde Parametrization to Original Parametrization --------------------
reparam_tilde_to_orig <- function(f_tilde, par_trans) {
  if (is.vector(f_tilde)) {
    f_orig <- f_tilde
    f_orig[par_trans == "logarithmic"] <- exp(f_tilde[par_trans == "logarithmic"])
    f_orig[par_trans == "logit"] <- 1 / (1 + exp(-f_tilde[par_trans == "logit"]))
  } else if (is.matrix(f_tilde)) {
    f_orig <- f_tilde
    f_orig[, par_trans == "logarithmic"] <- exp(f_tilde[, par_trans == "logarithmic"])
    f_orig[, par_trans == "logit"] <- 1 / (1 + exp(-f_tilde[, par_trans == "logit"]))
  }
  return(f_orig)
}
# ------------------------------------------------------------------------------


# Convert Tilde Parametrization to Original Parametrization as String ----------
reparam_tilde_to_orig_str <- function(par_trans) {
  str_vec <- c()
  for (i in 1:length(par_trans)) {
    if (par_trans[i] == "identity") {
      str_vec[i] <- paste0("f[, ", i ,"]")
    } else if (par_trans[i] == "logarithmic") {
      str_vec[i] <- paste0("exp(f[, ", i ,"])")
    } else if (par_trans[i] == "logit") {
      str_vec[i] <- paste0("1 / (1 + exp(-f[, ", i ,"]))")
    }
  }
  str_f_orig <- paste0("cbind(", paste(str_vec, collapse = ", "), ")")
  return(str_f_orig)
}
# ------------------------------------------------------------------------------


# Convert Original Parametrization to Tilde Parametrization --------------------
reparam_orig_to_tilde <- function(f_orig, par_trans) {
  if (is.vector(f_orig)) {
    f_tilde <- f_orig
    f_tilde[par_trans == "logarithmic"] <- log(f_orig[par_trans == "logarithmic"])
    f_tilde[par_trans == "logit"] <- log(f_orig[par_trans == "logit"] / (1 - f_orig[par_trans == "logit"]))
  } else if (is.matrix(f_orig)) {
    f_tilde <- f_orig
    f_tilde[, par_trans == "logarithmic"] <- log(f_orig[, par_trans == "logarithmic"])
    f_tilde[, par_trans == "logit"] <- log(f_orig[, par_trans == "logit"] / (1 - f_orig[, par_trans == "logit"]))
  }
  return(f_tilde)
}
# ------------------------------------------------------------------------------


# Compute Jacobian of Link Function --------------------------------------------
reparam_link_jacob <- function(f_orig, par_trans) {
  h_diag <- rep(1L, length(par_trans))
  h_diag[par_trans == "logarithmic"] <- 1 / f_orig[par_trans == "logarithmic"]
  h_diag[par_trans == "logit"] <- 1 / (f_orig[par_trans == "logit"] - f_orig[par_trans == "logit"]^2)
  h_matrix <- diag(h_diag, nrow = length(par_trans))
  return(h_matrix)
}
# ------------------------------------------------------------------------------


# Compute Jacobian of Link Function as String ----------------------------------
reparam_link_jacob_str <- function(par_trans) {
  str_vec <- c()
  for (i in 1:length(par_trans)) {
    if (par_trans[i] == "identity") {
      str_vec[i] <- paste0("1")
    } else if (par_trans[i] == "logarithmic") {
      str_vec[i] <- paste0("exp(-f[, ", i ,"])")
    } else if (par_trans[i] == "logit") {
      str_vec[i] <- paste0("(exp(f[, ", i ,"]) + exp(-f[, ", i ,"]) + 2)")
    }
  }
  str_h_matrix <- paste0("diag(c(", paste(str_vec, collapse = ", "), "), nrow = ", length(par_trans), ")")
  return(str_h_matrix)
}
# ------------------------------------------------------------------------------


# Compute Jacobian Inverse of Link Function ------------------------------------
reparam_link_jacob_inv <- function(f_orig, par_trans) {
  h_diag <- rep(1L, length(par_trans))
  h_diag[par_trans == "logarithmic"] <- f_orig[par_trans == "logarithmic"]
  h_diag[par_trans == "logit"] <- f_orig[par_trans == "logit"] - f_orig[par_trans == "logit"]^2
  h_matrix <- diag(h_diag, nrow = length(par_trans))
  return(h_matrix)
}
# ------------------------------------------------------------------------------


# Compute Jacobian Inverse of Link Function as String --------------------------
reparam_link_jacob_inv_str <- function(par_trans) {
  str_vec <- c()
  for (i in 1:length(par_trans)) {
    if (par_trans[i] == "identity") {
      str_vec[i] <- paste0("1")
    } else if (par_trans[i] == "logarithmic") {
      str_vec[i] <- paste0("exp(f[, ", i ,"])")
    } else if (par_trans[i] == "logit") {
      str_vec[i] <- paste0("exp(f[, ", i ,"]) / (1 + exp(f[, ", i ,"]))^2")
    }
  }
  h_matrix_str <- paste0("diag(c(", paste(str_vec, collapse = ", "), "), nrow = ", length(par_trans), ")")
  return(h_matrix_str)
}
# ------------------------------------------------------------------------------


