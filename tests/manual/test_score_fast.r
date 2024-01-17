
test_that("score_fast", {

  y_real <- -1.1
  y_duration <- 1.3
  y_interval <- 0.6
  y_count <- 3L

  par_mean <- -1.2
  par_scale <- 0.9
  par_shape1 <- 0.8
  par_shape2 <- 1.4
  par_df <- 4.3

  fun_fast <- setup_fun_score(distr = "alaplace", param = "meanscale", scaling = "unit", orthog = FALSE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "alaplace", param = "meanscale", scaling = "unit", orthog = FALSE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "alaplace", param = "meanscale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "alaplace", param = "meanscale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "alaplace", param = "meanscale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "alaplace", param = "meanscale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "unit", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(0), f = cbind(par_shape1))
  fun_slow <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "unit", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(0), f = cbind(par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(0), f = cbind(par_shape1))
  fun_slow <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(0), f = cbind(par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(0), f = cbind(par_shape1))
  fun_slow <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(0), f = cbind(par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "unit", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(1), f = cbind(par_shape1))
  fun_slow <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "unit", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(1), f = cbind(par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(1), f = cbind(par_shape1))
  fun_slow <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(1), f = cbind(par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(1), f = cbind(par_shape1))
  fun_slow <- setup_fun_score(distr = "bernoulli", param = "prob", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logit"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(1), f = cbind(par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "exp", param = "scale", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "exp", param = "scale", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "exp", param = "scale", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "exp", param = "scale", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "exp", param = "scale", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "exp", param = "scale", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "gamma", param = "scale", scaling = "unit", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "gamma", param = "scale", scaling = "unit", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "gamma", param = "scale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "gamma", param = "scale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "gamma", param = "scale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "gamma", param = "scale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "gengamma", param = "scale", scaling = "unit", orthog = FALSE, par_trans = c("logarithmic", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1, par_shape2))
  fun_slow <- setup_fun_score(distr = "gengamma", param = "scale", scaling = "unit", orthog = FALSE, par_trans = c("logarithmic", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1, par_shape2))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "gengamma", param = "scale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("logarithmic", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1, par_shape2))
  fun_slow <- setup_fun_score(distr = "gengamma", param = "scale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("logarithmic", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1, par_shape2))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "gengamma", param = "scale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("logarithmic", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1, par_shape2))
  fun_slow <- setup_fun_score(distr = "gengamma", param = "scale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("logarithmic", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1, par_shape2))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "geom", param = "mean", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "geom", param = "mean", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "geom", param = "mean", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "geom", param = "mean", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "geom", param = "mean", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "geom", param = "mean", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "laplace", param = "meanscale", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "laplace", param = "meanscale", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "laplace", param = "meanscale", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "laplace", param = "meanscale", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "laplace", param = "meanscale", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "laplace", param = "meanscale", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "logistic", param = "meanscale", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "logistic", param = "meanscale", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "logistic", param = "meanscale", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "logistic", param = "meanscale", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "logistic", param = "meanscale", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "logistic", param = "meanscale", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "negbin", param = "nb2", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "negbin", param = "nb2", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "negbin", param = "nb2", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "negbin", param = "nb2", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "negbin", param = "nb2", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "negbin", param = "nb2", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "logarithmic"), par_static = c(TRUE, FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "logarithmic"), par_static = c(TRUE, FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "logarithmic"), par_static = c(TRUE, FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "logarithmic"), par_static = c(TRUE, FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "logarithmic"), par_static = c(TRUE, FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "norm", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "logarithmic"), par_static = c(TRUE, FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "pois", param = "mean", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "pois", param = "mean", scaling = "unit", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "pois", param = "mean", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "pois", param = "mean", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "pois", param = "mean", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_count), f = cbind(par_scale))
  fun_slow <- setup_fun_score(distr = "pois", param = "mean", scaling = "fisher_inv", orthog = TRUE, par_trans = c("logarithmic"), par_static = c(FALSE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_count), f = cbind(par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "t", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  fun_slow <- setup_fun_score(distr = "t", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  fun_slow <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  fun_slow <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity", "identity"), par_static = c(FALSE, TRUE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "t", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "logarithmic", "identity"), par_static = c(TRUE, FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  fun_slow <- setup_fun_score(distr = "t", param = "meanvar", scaling = "unit", orthog = TRUE, par_trans = c("identity", "logarithmic", "identity"), par_static = c(TRUE, FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "logarithmic", "identity"), par_static = c(TRUE, FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  fun_slow <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "logarithmic", "identity"), par_static = c(TRUE, FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "logarithmic", "identity"), par_static = c(TRUE, FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  fun_slow <- setup_fun_score(distr = "t", param = "meanvar", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "logarithmic", "identity"), par_static = c(TRUE, FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale, par_df))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "vonmises", param = "meanconc", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "vonmises", param = "meanconc", scaling = "unit", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "vonmises", param = "meanconc", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "vonmises", param = "meanconc", scaling = "fisher_inv_sqrt", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "vonmises", param = "meanconc", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  fun_slow <- setup_fun_score(distr = "vonmises", param = "meanconc", scaling = "fisher_inv", orthog = TRUE, par_trans = c("identity", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_real), f = cbind(par_mean, par_scale))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "weibull", param = "scale", scaling = "unit", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "weibull", param = "scale", scaling = "unit", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "weibull", param = "scale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "weibull", param = "scale", scaling = "fisher_inv_sqrt", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

  fun_fast <- setup_fun_score(distr = "weibull", param = "scale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = TRUE)
  ss_fast <- fun_fast(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  fun_slow <- setup_fun_score(distr = "weibull", param = "scale", scaling = "fisher_inv", orthog = FALSE, par_trans = c("logarithmic", "identity"), par_static = c(FALSE, TRUE), fast = FALSE)
  ss_slow <- fun_slow(y = as.matrix(y_duration), f = cbind(par_scale, par_shape1))
  expect_equal(unname(ss_fast), unname(ss_slow), tolerance = 1e-9)

})
