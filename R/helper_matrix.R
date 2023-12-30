
# MATRIX FUNCTIONS


# Compute Matrix Inverse -------------------------------------------------------
matrix_inv <- function(mat) {
  mat <- as.matrix(mat)
  mat_inv <- try_and_be_silent(pracma::pinv(mat))
  if ("try-error" %in% class(mat_inv)) {
    mat_inv <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  }
  return(mat_inv)
}
# ------------------------------------------------------------------------------


# Compute Matrix Inverse Square Root -------------------------------------------
matrix_inv_sqrt <- function(mat) {
  mat <- as.matrix(mat)
  mat_eigen <- try_and_be_silent(base::eigen(mat, symmetric = TRUE))
  if ("try-error" %in% class(mat_eigen)) {
    mat_inv_sqrt <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  } else {
    mat_eigen$values[mat_eigen$values < 1e-6] <- 0
    mat_inv_sqrt <- try_and_be_silent(mat_eigen$vectors %*% pracma::pinv(mat_eigen$vectors %*% diag(sqrt(mat_eigen$values))))
    if ("try-error" %in% class(mat_inv_sqrt)) {
      mat_inv_sqrt <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
    }
  }
  return(mat_inv_sqrt)
}
# ------------------------------------------------------------------------------


# Compute Diagonal Matrix Inverse ----------------------------------------------
matrix_diag_inv <- function(mat) {
  mat <- as.matrix(mat)
  mat_inv <- be_silent(diag(1 / diag(mat), nrow = nrow(mat)))
  return(mat_inv)
}
# ------------------------------------------------------------------------------


# Compute Diagonal Matrix Inverse Square Root ----------------------------------
matrix_diag_inv_sqrt <- function(mat) {
  mat <- as.matrix(mat)
  mat_inv_sqrt <- be_silent(diag(1 / sqrt(diag(mat)), nrow = nrow(mat)))
  return(mat_inv_sqrt)
}
# ------------------------------------------------------------------------------


# Compute Partial Diagonal Matrix Inverse --------------------------------------
matrix_part_inv <- function(mat, idx) {
  mat <- as.matrix(mat)
  sub_mat_inv <- try_and_be_silent(pracma::pinv(mat[!idx, !idx]))
  if ("try-error" %in% class(sub_mat_inv)) {
    mat_inv <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  } else {
    mat_inv <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
    mat_inv[!idx, !idx] <- sub_mat_inv
    mat_inv[idx, idx] <- be_silent(diag(1 / diag(mat[idx, idx]), nrow = sum(idx)))
  }
  return(mat_inv)
}
# ------------------------------------------------------------------------------


# Compute Partial Diagonal Matrix Inverse Square Root --------------------------
matrix_part_inv_sqrt <- function(mat, idx) {
  mat <- as.matrix(mat)
  sub_mat_eigen <- try_and_be_silent(base::eigen(mat[!idx, !idx], symmetric = TRUE))
  if ("try-error" %in% class(sub_mat_eigen)) {
    mat_inv_sqrt <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  } else {
    sub_mat_eigen$values[sub_mat_eigen$values < 1e-6] <- 0
    sub_mat_inv_sqrt <- try_and_be_silent(sub_mat_eigen$vectors %*% pracma::pinv(sub_mat_eigen$vectors %*% diag(sqrt(sub_mat_eigen$values))))
    if ("try-error" %in% class(sub_mat_inv_sqrt)) {
      mat_inv_sqrt <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
    } else {
      mat_inv_sqrt <- matrix(0, nrow = nrow(mat), ncol = ncol(mat))
      mat_inv_sqrt[!idx, !idx] <- sub_mat_inv_sqrt
      mat_inv_sqrt[idx, idx] <- be_silent(diag(1 / sqrt(diag(mat[idx, idx])), nrow = sum(idx)))
    }
  }
  return(mat_inv_sqrt)
}
# ------------------------------------------------------------------------------


