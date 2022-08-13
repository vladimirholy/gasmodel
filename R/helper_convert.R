
# CONVERSION FUNCTIONS


# Convert Variance-Covariance Array to Variance Matrix -------------------------
convert_varcov_array_to_var_matrix <- function(varcov_array) {
  if (dim(varcov_array)[2] == 1L && dim(varcov_array)[3] == 1L) {
    var_matrix <- as.matrix(varcov_array[, 1, 1])
  } else {
    var_matrix <- t(apply(varcov_array, MARGIN = 1, FUN = diag))
  }
  return(var_matrix)
}
# ------------------------------------------------------------------------------


# Convert Variance-Covariance Matrix to Variance-Covariance Vector -------------
convert_varcov_matrix_to_varcov_vector <- function(varcov_mat) {
  n <- nrow(varcov_mat)
  varcov_vec <- rep(NA_real_, times = n + n * (n - 1) / 2)
  varcov_vec[1:n] <- diag(varcov_mat)
  varcov_vec[(n + 1):(n + n * (n - 1) / 2)] <- varcov_mat[lower.tri(varcov_mat)]
  return(varcov_vec)
}
# ------------------------------------------------------------------------------


# Convert Variance-Covariance Vector to Variance-Covariance Matrix -------------
convert_varcov_vector_to_varcov_matrix <- function(varcov_vec) {
  n <- sqrt(1 / 4 + 2 * length(varcov_vec)) - 1 / 2
  varcov_mat <- diag(varcov_vec[1:n])
  varcov_mat[lower.tri(varcov_mat)] <- varcov_vec[(n + 1):(n + n * (n - 1) / 2)]
  varcov_mat[upper.tri(varcov_mat)] <- t(varcov_mat)[upper.tri(varcov_mat)]
  return(varcov_mat)
}
# ------------------------------------------------------------------------------


# Convert Variance-Covariance Matrix to Parameter Vector -----------------------
convert_varcov_matrix_to_sc_vector <- function(varcov_mat) {
  n <- nrow(varcov_mat)
  sc_vec <- rep(NA_real_, times = n + n * (n - 1) / 2)
  sc_vec[1:n] <- diag(varcov_mat)
  sc_vec[(n + 1):(n + n * (n - 1) / 2)] <- varcov_mat[lower.tri(varcov_mat)] * 2
  return(sc_vec)
}
# ------------------------------------------------------------------------------


# Convert Parameter Vector to Variance-Covariance Matrix -----------------------
convert_sc_vector_to_varcov_matrix <- function(sc_vec) {
  n <- sqrt(1 / 4 + 2 * length(sc_vec)) - 1 / 2
  varcov_mat <- diag(sc_vec[1:n])
  varcov_mat[lower.tri(varcov_mat)] <- sc_vec[(n + 1):(n + n * (n - 1) / 2)] / 2
  varcov_mat[upper.tri(varcov_mat)] <- t(varcov_mat)[upper.tri(varcov_mat)]
  return(varcov_mat)
}
# ------------------------------------------------------------------------------


# Convert Kronecker Matrix to Parameter Matrix ---------------------------------
convert_krone_matrix_to_sc_matrix <- function(krone_matrix) {
  n <- sqrt(nrow(krone_matrix))
  sc_matrix <- matrix(nrow = n + n * (n - 1 ) / 2, ncol = n + n * (n - 1 ) / 2)
  sc_idx_var <- 1:n
  sc_idx_cov <- (n + 1):(n + n * (n - 1) / 2)
  krone_idx_var <- which(as.vector(as.logical(diag(n))))
  krone_idx_cov_lower <- matrix(1:n^2, nrow = n, ncol = n)[lower.tri(diag(n))]
  krone_idx_cov_upper <- t(matrix(1:n^2, nrow = n, ncol = n))[lower.tri(diag(n))]
  sc_matrix[sc_idx_var, sc_idx_var] <- krone_matrix[krone_idx_var, krone_idx_var]
  sc_matrix[sc_idx_var, sc_idx_cov] <- krone_matrix[krone_idx_var, krone_idx_cov_lower] + krone_matrix[krone_idx_var, krone_idx_cov_upper]
  sc_matrix[sc_idx_cov, sc_idx_var] <- krone_matrix[krone_idx_cov_lower, krone_idx_var] + krone_matrix[krone_idx_cov_upper, krone_idx_var]
  sc_matrix[sc_idx_cov, sc_idx_cov] <- krone_matrix[krone_idx_cov_lower, krone_idx_cov_lower] + krone_matrix[krone_idx_cov_lower, krone_idx_cov_upper] + krone_matrix[krone_idx_cov_upper, krone_idx_cov_lower] + krone_matrix[krone_idx_cov_upper, krone_idx_cov_upper]
  return(sc_matrix)
}
# ------------------------------------------------------------------------------


# Convert Parameter Matrix to Kronecker Matrix ---------------------------------
convert_sc_matrix_to_krone_matrix <- function(sc_matrix) {
  n <- sqrt(1 / 4 + 2 * nrow(sc_matrix)) - 1 / 2
  krone_matrix <- matrix(nrow = n * n, ncol = n * n)
  sc_idx_var <- 1:n
  sc_idx_cov <- (n + 1):(n + n * (n - 1) / 2)
  krone_idx_var <- which(as.vector(as.logical(diag(n))))
  krone_idx_cov_lower <- matrix(1:n^2, nrow = n, ncol = n)[lower.tri(diag(n))]
  krone_idx_cov_upper <- t(matrix(1:n^2, nrow = n, ncol = n))[lower.tri(diag(n))]
  krone_matrix[krone_idx_var, krone_idx_var] <- sc_matrix[sc_idx_var, sc_idx_var]
  krone_matrix[krone_idx_var, krone_idx_cov_lower] <- sc_matrix[sc_idx_var, sc_idx_cov] / 2
  krone_matrix[krone_idx_var, krone_idx_cov_upper] <- sc_matrix[sc_idx_var, sc_idx_cov] / 2
  krone_matrix[krone_idx_cov_lower, krone_idx_var] <- sc_matrix[sc_idx_cov, sc_idx_var] / 2
  krone_matrix[krone_idx_cov_upper, krone_idx_var] <- sc_matrix[sc_idx_cov, sc_idx_var] / 2
  krone_matrix[krone_idx_cov_lower, krone_idx_cov_lower] <- sc_matrix[sc_idx_cov, sc_idx_cov] / 4
  krone_matrix[krone_idx_cov_lower, krone_idx_cov_upper] <- sc_matrix[sc_idx_cov, sc_idx_cov] / 4
  krone_matrix[krone_idx_cov_upper, krone_idx_cov_lower] <- sc_matrix[sc_idx_cov, sc_idx_cov] / 4
  krone_matrix[krone_idx_cov_upper, krone_idx_cov_upper] <- sc_matrix[sc_idx_cov, sc_idx_cov] / 4
  return(krone_matrix)
}
# ------------------------------------------------------------------------------


# Convert Coefficient Vector to Optimization Vector  ---------------------------
convert_coef_vector_to_theta_vector <- function(coef_vec, coef_fix_value, coef_fix_other) {
  theta_vec <- coef_vec[is.na(coef_fix_value)]
  return(theta_vec)
}
# ------------------------------------------------------------------------------


# Convert Optimization Vector to Coefficient Vector ----------------------------
convert_theta_vector_to_coef_vector <- function(theta_vec, coef_fix_value, coef_fix_other) {
  coef_vec <- coef_fix_value + as.vector(coef_fix_other[, is.na(coef_fix_value), drop = FALSE] %*% theta_vec)
  coef_vec[is.na(coef_fix_value)] <- theta_vec
  return(coef_vec)
}


# Convert Coefficient Matrix to Optimization Matrix ----------------------------
convert_coef_matrix_to_theta_matrix <- function(coef_mat, coef_fix_value, coef_fix_other) {
  theta_mat <- coef_mat[is.na(coef_fix_value), is.na(coef_fix_value), drop = FALSE]
  return(theta_mat)
}
# ------------------------------------------------------------------------------


# Convert Optimization Matrix to Coefficient Matrix ----------------------------
convert_theta_matrix_to_coef_matrix <- function(theta_mat, coef_fix_value, coef_fix_other) {
  cnv <- coef_fix_other
  cnv[is.na(coef_fix_value), is.na(coef_fix_value)] <- diag(sum(is.na(coef_fix_value)))
  cnv <- cnv[, is.na(coef_fix_value), drop = FALSE]
  coef_mat <- cnv %*% theta_mat %*% t(cnv)
  return(coef_mat)
}
# ------------------------------------------------------------------------------


# Convert Coefficient Vector to Structured List --------------------------------
convert_coef_vector_to_struc_list <- function(coef_vec, m, p, q, par_names, par_of_coef_names) {
  dyn_names <- c("omega", "beta", "alpha", "phi")
  struc_list <- split(coef_vec, f = par_of_coef_names)
  struc_list <- lapply(1:length(par_names), function(i) { split(struc_list[[i]], f = factor(rep(dyn_names, c(1L, m[i], p[i], q[i])), levels = dyn_names)) })
  names(struc_list) <- par_names
  return(struc_list)
}
# ------------------------------------------------------------------------------


# Convert Structured List to Coefficient Vector --------------------------------
convert_struc_list_to_coef_vector <- function(struc_list) {
  coef_vec <- unlist(struc_list)
  return(coef_vec)
}
# ------------------------------------------------------------------------------


# Convert Optimization Vector to Structured List -------------------------------
convert_theta_vector_to_struc_list <- function(theta_vec, coef_fix_value, coef_fix_other, m, p, q, par_names, par_of_coef_names) {
  coef_vec <- convert_theta_vector_to_coef_vector(theta_vec, coef_fix_value = coef_fix_value, coef_fix_other = coef_fix_other)
  struc_list <- convert_coef_vector_to_struc_list(coef_vec, m = m, p = p, q = q, par_names = par_names, par_of_coef_names = par_of_coef_names)
  return(struc_list)
}
# ------------------------------------------------------------------------------


# Convert Structured List to Optimization Vector -------------------------------
convert_struc_list_to_theta_vector <- function(struc_list, coef_fix_value, coef_fix_other) {
  coef_vec <- convert_struc_list_to_coef_vector(struc_list)
  theta_vec <- convert_coef_vector_to_theta_vector(coef_vec, coef_fix_value = coef_fix_value, coef_fix_other = coef_fix_other)
  return(theta_vec)
}
# ------------------------------------------------------------------------------


