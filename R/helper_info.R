
# MODEL INFORMATION FUNCTIONS


# Get the Index Names ----------------------------------------------------------
info_data <- function(y, x, skip_t = 0L) {
  info <- list()
  if (!is.null(rownames(y))) {
    info$index_time <- rownames(y)
  } else if (any(sapply(x, function(e) { !is.null(rownames(e)) }))) {
    info$index_time <- rownames(x[[which(sapply(x, function(e) { !is.null(rownames(e)) }))[1]]])
  } else {
    info$index_time <- paste0("t", (skip_t + 1):(skip_t + nrow(y)))
  }
  if (!is.null(colnames(y))) {
    info$index_series <- colnames(y)
  } else {
    info$index_series <- paste0("ser", 1:ncol(y))
  }
  info$index_time_list <- rep(list(info$index_time), times = length(x))
  info$index_vars_list <- list()
  for (i in 1:length(x)) {
    if (ncol(x[[i]]) == 0L) {
      info$index_vars_list[[i]] <- character()
    } else if (!is.null(colnames(x[[i]]))) {
      info$index_vars_list[[i]] <- colnames(x[[i]])
    } else {
      info$index_vars_list[[i]] <- paste0("var", 1:ncol(x[[i]]))
    }
  }
  return(info)
}
# ------------------------------------------------------------------------------


# Get the Distribution Characteristics -----------------------------------------
info_distribution <- function(distr, param) {
  details_distr <- distr(filter_distr = distr, filter_param = param)[1, ]
  info <- list()
  info$type <- details_distr$type
  info$dim <- details_distr$dim
  info$orthog <- details_distr$orthog
  info$default <- details_distr$default
  return(info)
}
# ------------------------------------------------------------------------------


# Get the Names and Numbers of Groups and Parameters ---------------------------
info_parameters <- function(distr, param, n) {
  details_par <- do.call(paste("distr", distr, param, "parameters", sep = "_"), args = list(n = n))
  info <- list()
  info$group_names <- unique(details_par$group_of_par_names)
  info$group_num <- length(info$group_names)
  info$group_of_par_names <- factor(details_par$group_of_par_names, levels = info$group_names)
  info$par_names <- details_par$par_names
  info$par_num <- length(info$par_names)
  info$par_in_group_num <- sapply(split(info$par_names, f = factor(info$group_of_par_names, levels = info$group_names)), length)
  info$par_support <- details_par$par_support
  info$par_trans <- rep("identity", times = info$par_num)
  return(info)
}
# ------------------------------------------------------------------------------


# Update the Logarithmic and Logistic Transformations --------------------------
info_linked_parameters <- function(info_par, par_link) {
  info <- info_par
  info$par_trans[par_link & info$par_support == "positive"] <- "logarithmic"
  info$par_trans[par_link & info$par_support == "probability"] <- "logit"
  info$par_support[par_link] <- "real"
  info$par_names[info$par_trans == "logarithmic"] <- paste0("log(", info$par_names[info$par_trans == "logarithmic"], ")")
  info$par_names[info$par_trans == "logit"] <- paste0("logit(", info$par_names[info$par_trans == "logit"], ")")
  return(info)
}
# ------------------------------------------------------------------------------


# Get the Names and the Number of Coefficients ---------------------------------
info_coefficients <- function(m, p, q, par_static, par_names, par_num, group_names, group_of_par_names) {
  info <- list()
  info$coef_names <- unlist(lapply(1:par_num, function(i) { if (par_static[i]) { par_names[i] } else { paste0(par_names[i], c("_omega", if (m[i] > 0L) { paste0("_beta", 1:m[i]) }, if (p[i] > 0L) { paste0("_alpha", 1:p[i]) }, if (q[i] > 0L) { paste0("_phi", 1:q[i]) })) } }))
  info$coef_num <- length(info$coef_names)
  info$coef_in_par_num <- 1L + m + p + q
  info$coef_in_group_num <- sapply(split(info$coef_in_par_num, f = factor(group_of_par_names, levels = group_names)), sum)
  info$par_of_coef_names <- factor(rep(par_names, times = info$coef_in_par_num), levels = par_names)
  info$dyn_names <- c("omega", "beta", "alpha", "phi")
  return(info)
}
# ------------------------------------------------------------------------------


# Get the Names and the Number of Thetas ---------------------------------------
info_thetas <- function(coef_fix_value, coef_fix_other, coef_names) {
  info <- list()
  info$theta_names <- convert_coef_vector_to_theta_vector(coef_names, coef_fix_value = coef_fix_value, coef_fix_other = coef_fix_other)
  info$theta_num <- length(info$theta_names)
  return(info)
}
# ------------------------------------------------------------------------------


# Get the Title of the GAS Model -----------------------------------------------
info_title <- function(distr, param, scaling) {
  info <- list()
  info$title_distr <- distr(filter_distr = distr, filter_param = param)$title[1]
  info$title_scaling <- switch(scaling, "unit" = "Unit Scaling", "fisher_inv" = "Fisher Inverse Scaling", "fisher_inv_sqrt" = "Fisher Square Root Inverse Scaling")
  info$title <- paste(info$title_distr, info$title_scaling, sep = " / ")
  return(info)
}
# ------------------------------------------------------------------------------


