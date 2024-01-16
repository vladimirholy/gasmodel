
test_that("distr_lognorm_logmeanvar", {

  f <- c(0.8, 1.3)
  by <- 1e-4
  y <- as.matrix(seq(from = by, to = 1e4, by = by))
  t <- nrow(y)
  f_all <- matrix(f, nrow = t, ncol = length(f), byrow = TRUE)

  density <- distr_lognorm_logmeanvar_density(y, f_all)
  expect_equal(sum(density * by), 1, tolerance = 1e-6)

  loglik <- distr_lognorm_logmeanvar_loglik(y, f_all)
  expect_equal(density, exp(loglik), tolerance = 1e-12)

  mean_num <- sum(density * by * y)
  mean <- as.vector(distr_lognorm_logmeanvar_mean(t(f)))
  expect_equal(mean, mean_num, tolerance = 1e-6)

  var_num <- sum(density * by * (y - mean_num)^2)
  var <- as.vector(distr_lognorm_logmeanvar_var(t(f)))
  expect_equal(var, var_num, tolerance = 1e-4)

  fisher_num <- matrix(0, nrow = length(f), ncol = length(f))
  for (i in 1:t) {
    score <- distr_lognorm_logmeanvar_score(y[i, , drop = FALSE], t(f))
    fisher_num <- fisher_num + t(score) %*% score * density[i, ] * by
  }
  fisher <- distr_lognorm_logmeanvar_fisher(t(f))[1, , ]
  expect_equal(fisher, fisher_num, tolerance = 1e-6)

  f_start <- distr_lognorm_logmeanvar_start(distr_lognorm_logmeanvar_random(t = 1e6, f = f))
  expect_equal(f, f_start, tolerance = 1e-2)

})
