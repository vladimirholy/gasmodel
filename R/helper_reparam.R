
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


# Compute Jacobian Inverse of Link Function ------------------------------------
reparam_link_jacob_inv <- function(f_orig, par_trans) {
  h_diag <- rep(1L, length(par_trans))
  h_diag[par_trans == "logarithmic"] <- f_orig[par_trans == "logarithmic"]
  h_diag[par_trans == "logit"] <- f_orig[par_trans == "logit"] - f_orig[par_trans == "logit"]^2
  h_matrix <- diag(h_diag, nrow = length(par_trans))
  return(h_matrix)
}
# ------------------------------------------------------------------------------


