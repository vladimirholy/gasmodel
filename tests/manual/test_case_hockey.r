
test_that("test_case_hockey", {

  data("ice_hockey_championships")

  y <- ice_hockey_championships$rankings[c(1:22, 24:25), ]

  est_gas <- gas(y, distr = "pluce", scaling = "unit", reg = "sep", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(est_gas, "gas")
  expect_equal(coef(est_gas)[1:9], setNames(c(-0.805, 0.44, 0.485, 0.227, 0.44, 0.485, 3.605, 0.44, 0.485), c("log(worth1)_omega", "log(worth1)_alpha1", "log(worth1)_phi1", "log(worth2)_omega", "log(worth2)_alpha1", "log(worth2)_phi1", "log(worth3)_omega", "log(worth3)_alpha1", "log(worth3)_phi1")), tolerance = 1e-3)

})
