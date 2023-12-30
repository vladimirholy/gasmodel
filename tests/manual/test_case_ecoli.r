
test_that("test_case_ecoli", {

  data("ecoli", package = "tscount")

  y <- ecoli$cases

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "unit", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.956, 0.037, 0.889, 0.057), c("log(mean)_omega", "log(mean)_alpha1", "log(mean)_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv_sqrt", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.953, 0.111, 0.894, 0.057), c("log(mean)_omega", "log(mean)_alpha1", "log(mean)_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.948, 0.335, 0.901, 0.057), c("log(mean)_omega", "log(mean)_alpha1", "log(mean)_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.467, 11.909, 0.987, 0.066), c("mean_omega", "mean_alpha1", "mean_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(19.872, 2.135, 0.886, 0.060), c("mean_omega", "mean_alpha1", "mean_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(19.828, 0.339, 0.874, 0.060), c("mean_omega", "mean_alpha1", "mean_phi1", "dispersion")), tolerance = 1e-3)

})
