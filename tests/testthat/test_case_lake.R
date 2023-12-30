
test_that("case_lake", {

  data("LakeHuron")

  y <- LakeHuron - 570
  x <- 1:length(y)

  est_gas <- gas(y = y, x = x, distr = "norm", regress = "sep", par_static = c(FALSE, TRUE))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas), setNames(c(9.989, -0.020, 0.464, 0.665, 0.457), c("mean_omega", "mean_beta1", "mean_alpha1", "mean_phi1", "var")), tolerance = 1e-3)

})
