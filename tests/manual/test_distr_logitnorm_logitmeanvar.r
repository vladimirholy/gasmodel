
test_that("distr_logitnorm_logitmeanvar", {

  f <- c(0.8, 1.3)
  by <- 1e-3
  y <- as.matrix(1 / (1 + exp(-seq(from = -1e2, to = 1e2, by = by))))
  t <- nrow(y)
  f_all <- matrix(f, nrow = t, ncol = length(f), byrow = TRUE)

  density <- distr_logitnorm_logitmeanvar_density(y, f_all)
  expect_equal(sum(density * by), 1, tolerance = 1e-6)

  loglik <- distr_logitnorm_logitmeanvar_loglik(y, f_all)
  expect_equal(density, exp(loglik), tolerance = 1e-12)

  mean_num <- sum(density * by * y)
  mean <- as.vector(distr_logitnorm_logitmeanvar_mean(t(f)))
  expect_equal(mean, mean_num, tolerance = 1e-3)

  var_num <- sum(density * by * (y - mean_num)^2)
  var <- as.vector(distr_logitnorm_logitmeanvar_var(t(f)))
  expect_equal(var, var_num, tolerance = 1e-2)

  fisher_num <- matrix(0, nrow = length(f), ncol = length(f))
  for (i in 1:t) {
    score <- distr_logitnorm_logitmeanvar_score(y[i, , drop = FALSE], t(f))
    if (all(is.finite(score))) {
      fisher_num <- fisher_num + t(score) %*% score * density[i, ] * by
    }
  }
  fisher <- distr_logitnorm_logitmeanvar_fisher(t(f))[1, , ]
  expect_equal(fisher, fisher_num, tolerance = 1e-6)

  f_start <- distr_logitnorm_logitmeanvar_start(distr_logitnorm_logitmeanvar_random(t = 1e6, f = f))
  expect_equal(f, f_start, tolerance = 1e-2)

})
