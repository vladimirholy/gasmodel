
# SPECIAL TYPES OF FIXED COEFFICIENTS FUNCTION


# Update Fixed Coefficients Using Special Structure ----------------------------
fixed_coefficients <- function(model, info_par, info_coef) {
  coef_fix_value <- model$coef_fix_value
  coef_fix_other <- model$coef_fix_other
  if ("panel_structure" %in% model$coef_fix_special) {
    for (i in 1:info_par$group_num) {
      if (length(unique(model$m[info_par$group_of_par_names == info_par$group_names[i]])) == 1L && length(unique(model$p[info_par$group_of_par_names == info_par$group_names[i]])) == 1L && length(unique(model$q[info_par$group_of_par_names == info_par$group_names[i]])) == 1L) {
        first_par <- which(info_par$group_of_par_names == info_par$group_names[i])[1]
        other_par <- which(info_par$group_of_par_names == info_par$group_names[i])[-1]
        first_par_dyn_coef <- which(info_coef$par_of_coef_names == info_par$par_names[first_par])[-1]
        for (j in 1:length(other_par)) {
          cur_par_dyn_coef <- which(info_coef$par_of_coef_names == info_par$par_names[other_par[j]])[-1]
          for (k in 1:length(cur_par_dyn_coef)) {
            coef_fix_value[cur_par_dyn_coef[k]] <- 0
            coef_fix_other[cur_par_dyn_coef[k], ] <- 0
            coef_fix_other[cur_par_dyn_coef[k], first_par_dyn_coef[k]] <- 1
          }
        }
      }
    }
  }
  if ("zero_sum_intercept" %in% model$coef_fix_special) {
    for (i in 1:info_par$group_num) {
      first_par_omega_coef <- sapply(which(info_par$group_of_par_names == info_par$group_names[i])[1], function(j) { which(info_coef$par_of_coef_names == info_par$par_names[j])[1] })
      other_par_omega_coef <- sapply(which(info_par$group_of_par_names == info_par$group_names[i])[-1], function(j) { which(info_coef$par_of_coef_names == info_par$par_names[j])[1] })
      coef_fix_value[first_par_omega_coef] <- 0
      coef_fix_other[first_par_omega_coef, ] <- 0
      coef_fix_other[first_par_omega_coef, other_par_omega_coef] <- -1
    }
  }
  if ("random_walk" %in% model$coef_fix_special) {
    for (j in 1:info_par$par_num) {
      if (model$q[j] > 0L) {
        cur_par_phi_coef <- which(info_coef$par_of_coef_names == info_par$par_names[j])[(1L + model$m[j] + model$p[j] + 1L):(1L + model$m[j] + model$p[j] + model$q[j])]
        coef_fix_value[cur_par_phi_coef] <- 1
        coef_fix_other[cur_par_phi_coef, ] <- 0
      }
    }
  }
  coef_fix_other[is.na(coef_fix_value), ] <- NA_real_
  coef_fix_other[, !is.na(coef_fix_value)] <- NA_real_
  coef_fix_value <- name_vector(coef_fix_value, info_coef$coef_names)
  coef_fix_other <- name_matrix(coef_fix_other, info_coef$coef_names, info_coef$coef_names)
  report <- list(coef_fix_value = coef_fix_value, coef_fix_other = coef_fix_other)
  return(report)
}


