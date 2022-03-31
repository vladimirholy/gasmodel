
# NAMING FUNCTIONS


# Set Vector Names -------------------------------------------------------------
name_vector <- function(obj, elem_names) {
  names(obj) <- elem_names
  return(obj)
}
# ------------------------------------------------------------------------------


# Set Matrix Names -------------------------------------------------------------
name_matrix <- function(obj, row_names, col_names, drop = c(FALSE, FALSE), zero = c(FALSE, FALSE)) {
  rownames(obj) <- row_names
  colnames(obj) <- col_names
  if (any(zero & dim(obj) == 0L)) {
    obj <- numeric()
  } else if (any(drop)) {
    obj <- abind::adrop(obj, drop = (drop & dim(obj) == 1L))
  }
  return(obj)
}
# ------------------------------------------------------------------------------


# Set Three-Dimensional Array Names --------------------------------------------
name_3d_array <- function(obj, d1_names, d2_names, d3_names, drop = c(FALSE, FALSE, FALSE), zero = c(FALSE, FALSE, FALSE)) {
  dimnames(obj) <- list(d1_names, d2_names, d3_names)
  if (any(zero & dim(obj) == 0L)) {
    obj <- numeric()
  } else if (any(drop)) {
    obj <- abind::adrop(obj, drop = (drop & dim(obj) == 1L))
  }
  return(obj)
}
# ------------------------------------------------------------------------------


# Set List of Matrices Names ---------------------------------------------------
name_list_of_matrices <- function(obj, elem_names, row_names_list, col_names_list, drop = c(FALSE, FALSE), zero = c(FALSE, FALSE)) {
  names(obj) <- elem_names
  for (i in 1:length(obj)) {
    rownames(obj[[i]]) <- row_names_list[[i]]
    colnames(obj[[i]]) <- col_names_list[[i]]
    if (any(zero & dim(obj[[i]]) == 0L)) {
      obj[[i]] <- numeric()
    } else if (any(drop)) {
      obj[[i]] <- abind::adrop(obj[[i]], drop = (drop & dim(obj[[i]]) == 1L))
    }
  }
  return(obj)
}
# ------------------------------------------------------------------------------


