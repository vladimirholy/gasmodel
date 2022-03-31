
# MATRIX FUNCTIONS


# Compute Matrix Inverse -------------------------------------------------------
matrix_inv <- function(mat) {
  mat_inv <- suppressWarnings(try(pracma::pinv(mat), silent = TRUE))
  if (class(mat_inv) == "try-error") {
    mat_inv <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  }
  return(mat_inv)
}
# ------------------------------------------------------------------------------


# Compute Matrix Inverse Square Root -------------------------------------------
matrix_inv_sqrt <- function(mat) {
  mat_eigen <- base::eigen(mat, symmetric = TRUE)
  mat_eigen$values[mat_eigen$values < 1e-6] <- 0
  mat_inv_sqrt <- suppressWarnings(try(mat_eigen$vectors %*% pracma::pinv(mat_eigen$vectors %*% diag(sqrt(mat_eigen$values))), silent = TRUE))
  if (class(mat_inv_sqrt) == "try-error") {
    mat_inv_sqrt <- matrix(NA_real_, nrow = nrow(mat), ncol = ncol(mat))
  }
  return(mat_inv_sqrt)
}
# ------------------------------------------------------------------------------


# Compute Diagonal Matrix Inverse ----------------------------------------------
matrix_diag_inv <- function(mat) {
  mat_inv <- suppressWarnings(diag(1 / diag(mat), nrow = nrow(mat)))
  return(mat_inv)
}
# ------------------------------------------------------------------------------


# Compute Diagonal Matrix Inverse Square Root ----------------------------------
matrix_diag_inv_sqrt <- function(mat) {
  mat_inv_sqrt <- suppressWarnings(diag(1 / sqrt(diag(mat)), nrow = nrow(mat)))
  return(mat_inv_sqrt)
}
# ------------------------------------------------------------------------------


