
test_that("hockey", {
  y <- ice_hockey_championships$rankings[1:22, ]
  x <- lapply(1:ncol(ice_hockey_championships$host), function(i) { ice_hockey_championships$host[1:22, i] })
  res_estimate <- gas(y = y, x = x, distr = "pluce", coef_fix_special = c("zero_sum_intercept", "panel_structure"))
  expect_s3_class(res_estimate, "gas")
  expect_equal(res_estimate$fit$coef_est[4], setNames(0.506, "log(worth1)_phi1"), tolerance = 1e-3)
})
