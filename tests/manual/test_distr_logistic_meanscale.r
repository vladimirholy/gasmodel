
test_that("distr_logistic_meanscale", {

  f <- c(-1.2, 2.3)
  by <- 1e-4
  y <- as.matrix(seq(from = -1e3, to = 1e3, by = by))
  t <- nrow(y)
  f_all <- matrix(f, nrow = t, ncol = length(f), byrow = TRUE)

  density <- distr_logistic_meanscale_density(y, f_all)
  expect_equal(sum(density * by), 1, tolerance = 1e-6)

  loglik <- distr_logistic_meanscale_loglik(y, f_all)
  expect_equal(density, exp(loglik), tolerance = 1e-12)

  mean_num <- sum(density * by * y)
  mean <- as.vector(distr_logistic_meanscale_mean(t(f)))
  expect_equal(mean, mean_num, tolerance = 1e-6)

  var_num <- sum(density * by * (y - mean_num)^2)
  var <- as.vector(distr_logistic_meanscale_var(t(f)))
  expect_equal(var, var_num, tolerance = 1e-6)

  fisher_num <- matrix(0, nrow = length(f), ncol = length(f))
  for (i in 1:t) {
    score <- distr_logistic_meanscale_score(y[i, , drop = FALSE], t(f))
    fisher_num <- fisher_num + t(score) %*% score * density[i, ] * by
  }
  fisher <- distr_logistic_meanscale_fisher(t(f))[1, , ]
  expect_equal(fisher, fisher_num, tolerance = 1e-3)

  f_start <- distr_logistic_meanscale_start(distr_logistic_meanscale_random(t = 1e6, f = f))
  expect_equal(f, f_start, tolerance = 1e-1)

})
