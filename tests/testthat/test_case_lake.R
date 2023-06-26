
test_that("case_lake", {

  data(LakeHuron)

  y <- LakeHuron - 570
  x <- 1:length(y)

  est_gas <- gas(y = y, x = x, distr = "norm", regress = "sep", par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(est_gas$fit$coef_est[4], setNames(0.665, "mean_phi1"), tolerance = 1e-3)

})
