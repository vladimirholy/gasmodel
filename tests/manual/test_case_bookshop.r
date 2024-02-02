
test_that("test_case_bookshop", {

  data("bookshop_orders", package = "gasmodel")

  y <- bookshop_orders$duration_adj[-1]

  est_gas <- gas(y, distr = "exp", reg = "sep", scaling = "unit")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.024, 0.050, 0.963), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "exp", reg = "sep", scaling = "fisher_inv_sqrt")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.024, 0.050, 0.963), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "exp", reg = "sep", scaling = "fisher_inv")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.024, 0.050, 0.963), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gamma", reg = "sep", scaling = "unit")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.028, 0.053, 0.963, 0.949), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gamma", reg = "sep", scaling = "fisher_inv_sqrt")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.028, 0.051, 0.963, 0.949), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gamma", reg = "sep", scaling = "fisher_inv")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.028, 0.050, 0.963, 0.949), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gengamma", reg = "sep", scaling = "unit")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-1.161, 0.072, 0.950, 1.886, 0.661), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape1", "shape2")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gengamma", reg = "sep", scaling = "fisher_inv_sqrt")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-1.161, 0.065, 0.950, 1.886, 0.661), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape1", "shape2")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gengamma", reg = "sep", scaling = "fisher_inv")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-1.161, 0.059, 0.950, 1.886, 0.661), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape1", "shape2")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "unit", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.050, 0.057, 0.962, 0.947), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.049, 0.061, 0.957, 0.949), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.049, 0.057, 0.956, 0.948), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.967, 0.051, 0.953, 0.947), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "diag_fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.972, 0.058, 0.954, 0.948), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "diag_fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.977, 0.058, 0.960, 0.947), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.972, 0.064, 0.947, 0.948), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.975, 0.058, 0.955, 0.948), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

})
