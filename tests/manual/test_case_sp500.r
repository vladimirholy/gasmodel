
test_that("test_case_sp500", {

  data("sp500ret", package = "rugarch")

  y <- sp500ret$SP500RET

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "unit", par_link = c(FALSE, TRUE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.001, -9.587, 0.157, 0.993, 6.606), c("mean", "log(var)_omega", "log(var)_alpha1", "log(var)_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, TRUE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.001, -9.562, 0.092, 0.993, 6.778), c("mean", "log(var)_omega", "log(var)_alpha1", "log(var)_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, TRUE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, -8.802, 0.063, 0.995, 5.731), c("mean", "log(var)_omega", "log(var)_alpha1", "log(var)_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, 0, 0, 0.999, 4.104), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "diag_fisher_inv_sqrt", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.001, 0, 0, 0.993, 5.715), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "diag_fisher_inv", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.001, 0, 0.056, 0.997, 6.567), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.001, 0, 0, 0.993, 5.715), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.001, 0, 0.043, 1, 10.599), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

})
