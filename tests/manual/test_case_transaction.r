
test_that("test_case_transaction", {

  data("adjDurData", package = "ACDm")

  y <- adjDurData %>%
    dplyr::filter(as.Date(time) == "2009-05-15") %>%
    dplyr::select(adjDur) %>%
    as.matrix() %>%
    as.vector()

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "unit", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.014, 0.059, 0.993, 0.921), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.041, 0.055, 0.995, 0.92), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.063, 0.049, 0.995, 0.915), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(1.072, 0.072, 0.895, 0.893), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "diag_fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(1.1, 0.05, 0.989, 0.922), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "diag_fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(1.107, 0.049, 0.997, 0.92), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(1.096, 0.05, 0.991, 0.921), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(1.023, 0.046, 0.998, 0.911), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

})
