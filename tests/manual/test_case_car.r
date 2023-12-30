
test_that("test_case_car", {

  data("german_car_market_cap")

  y <- german_car_market_cap %>%
    dplyr::group_by(year) %>%
    dplyr::mutate(rank = rank(-market_cap)) %>%
    dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols = year, names_from = car_manufacturer, values_from = rank) %>%
    tibble::column_to_rownames("year") %>%
    as.matrix()

  est_gas <- gas(y, distr = "pluce", scaling = "unit", reg = "sep", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.571, 1.291, 0.695, 0.426, 1.291, 0.695, -0.997, 1.291, 0.695), c("log(worth1)_omega", "log(worth1)_alpha1", "log(worth1)_phi1", "log(worth2)_omega", "log(worth2)_alpha1", "log(worth2)_phi1", "log(worth3)_omega", "log(worth3)_alpha1", "log(worth3)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "pluce", scaling = "full_fisher_inv_sqrt", reg = "sep", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.697, 0.865, 0.696, 0.257, 0.865, 0.696, -0.953, 0.865, 0.696), c("log(worth1)_omega", "log(worth1)_alpha1", "log(worth1)_phi1", "log(worth2)_omega", "log(worth2)_alpha1", "log(worth2)_phi1", "log(worth3)_omega", "log(worth3)_alpha1", "log(worth3)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "pluce", scaling = "full_fisher_inv", reg = "sep", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.106, 0.584, 0.863, -1.661, 0.584, 0.863, -0.444, 0.584, 0.863), c("log(worth1)_omega", "log(worth1)_alpha1", "log(worth1)_phi1", "log(worth2)_omega", "log(worth2)_alpha1", "log(worth2)_phi1", "log(worth3)_omega", "log(worth3)_alpha1", "log(worth3)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "pluce", scaling = "diag_fisher_inv_sqrt", reg = "sep", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.617, 0.671, 0.699, 0.297, 0.671, 0.699, -0.914, 0.671, 0.699), c("log(worth1)_omega", "log(worth1)_alpha1", "log(worth1)_phi1", "log(worth2)_omega", "log(worth2)_alpha1", "log(worth2)_phi1", "log(worth3)_omega", "log(worth3)_alpha1", "log(worth3)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "pluce", scaling = "diag_fisher_inv", reg = "sep", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.573, 0.357, 0.691, 0.321, 0.357, 0.691, -0.894, 0.357, 0.691), c("log(worth1)_omega", "log(worth1)_alpha1", "log(worth1)_phi1", "log(worth2)_omega", "log(worth2)_alpha1", "log(worth2)_phi1", "log(worth3)_omega", "log(worth3)_alpha1", "log(worth3)_phi1")), tolerance = 1e-3)

})
