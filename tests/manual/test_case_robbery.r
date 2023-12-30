
test_that("test_case_robbery", {

  data("RobberyConvict", package = "glarma")

  y <- RobberyConvict[, 9]

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "unit", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.638, 0.025, 0.980, 0.043), c("log(mean)_omega", "log(mean)_alpha1", "log(mean)_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv_sqrt", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.628, 0.080, 0.981, 0.042), c("log(mean)_omega", "log(mean)_alpha1", "log(mean)_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(2.612, 0.257, 0.983, 0.043), c("log(mean)_omega", "log(mean)_alpha1", "log(mean)_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(13.940, 6.976, 0.979, 0.043), c("mean_omega", "mean_alpha1", "mean_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(14.357, 1.402, 0.977, 0.044), c("mean_omega", "mean_alpha1", "mean_phi1", "dispersion")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "negbin", reg = "sep", scaling = "fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(14.570, 0.246, 0.976, 0.043), c("mean_omega", "mean_alpha1", "mean_phi1", "dispersion")), tolerance = 1e-3)

})
