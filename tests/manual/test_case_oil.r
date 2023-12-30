
test_that("test_case_oil", {

  data("oil", package = "gamlss.data")

  y <- diff(oil$OILPRICE)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "unit", par_link = c(FALSE, TRUE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, -8.58, 0.188, 0.993, 7.421), c("mean", "log(var)_omega", "log(var)_alpha1", "log(var)_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, TRUE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, -8.567, 0.113, 0.992, 7.498), c("mean", "log(var)_omega", "log(var)_alpha1", "log(var)_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, TRUE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, -8.356, 0.059, 0.992, 12.842), c("mean", "log(var)_omega", "log(var)_alpha1", "log(var)_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.002, 0, 0, 1, 9.531), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "diag_fisher_inv_sqrt", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, 0, 0, 1, 6.471), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "diag_fisher_inv", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, 0, 0.071, 0.992, 7.393), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, 0, 0, 0.995, 6.654), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "t", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, FALSE, FALSE), par_static = c(TRUE, FALSE, TRUE), coef_bound_upper = c(Inf, Inf, Inf, 0.9999, Inf))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0, 0, 0.048, 1, 8), c("mean", "var_omega", "var_alpha1", "var_phi1", "df")), tolerance = 1e-3)

})
