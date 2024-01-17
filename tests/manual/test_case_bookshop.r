
test_that("test_case_bookshop", {

  data("bookshop_sales", package = "gasmodel")

  data_dur <- bookshop_sales %>%
    dplyr::mutate(duration = as.numeric(time - dplyr::lag(time)) / 60) %>%
    dplyr::mutate(duration = dplyr::recode(duration, "0" = 0.5)) %>%
    dplyr::mutate(just_time = as.vector(as.POSIXct(format(time, "1970-01-01 %H:%M:%S"), tz = "UTC"))) %>%
    tidyr::drop_na()

  model_spline <- stats::smooth.spline(data_dur$just_time, data_dur$duration, df = 10)

  data_dur <- data_dur %>%
    dplyr::mutate(duration_spline = stats::predict(model_spline, x = just_time)$y) %>%
    dplyr::mutate(duration_adj = duration / duration_spline) %>%
    dplyr::select(-just_time)

  y <- data_dur$duration_adj

  est_gas <- gas(y, distr = "exp", reg = "sep", scaling = "unit")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.023, 0.049, 0.963), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "exp", reg = "sep", scaling = "fisher_inv_sqrt")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.023, 0.049, 0.963), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "exp", reg = "sep", scaling = "fisher_inv")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.023, 0.049, 0.963), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gamma", reg = "sep", scaling = "unit")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.036, 0.052, 0.963, 0.942), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gamma", reg = "sep", scaling = "fisher_inv_sqrt")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.036, 0.050, 0.963, 0.942), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gamma", reg = "sep", scaling = "fisher_inv")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.036, 0.049, 0.963, 0.942), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gengamma", reg = "sep", scaling = "unit")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-1.019, 0.070, 0.952, 1.764, 0.683), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape1", "shape2")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gengamma", reg = "sep", scaling = "fisher_inv_sqrt")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-1.019, 0.063, 0.952, 1.764, 0.683), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape1", "shape2")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "gengamma", reg = "sep", scaling = "fisher_inv")
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-1.019, 0.057, 0.952, 1.764, 0.683), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape1", "shape2")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "unit", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.051, 0.056, 0.962, 0.944), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.05, 0.059, 0.959, 0.946), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv", par_link = c(TRUE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(-0.049, 0.054, 0.959, 0.945), c("log(scale)_omega", "log(scale)_alpha1", "log(scale)_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "unit", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.966, 0.049, 0.955, 0.944), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "diag_fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.971, 0.057, 0.954, 0.945), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "diag_fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.975, 0.057, 0.96, 0.945), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv_sqrt", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.971, 0.062, 0.949, 0.945), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

  est_gas <- gas(y, distr = "weibull", reg = "sep", scaling = "full_fisher_inv", par_link = c(FALSE, FALSE), par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(0.975, 0.055, 0.959, 0.945), c("scale_omega", "scale_alpha1", "scale_phi1", "shape")), tolerance = 1e-3)

})
