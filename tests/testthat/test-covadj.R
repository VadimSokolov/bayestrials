test_that("ANCOVA gain is 1 - r2", {
  g <- covadj_ancova_gain(0.3, n_unadjusted = 100)
  expect_equal(g$var_factor, 0.7)
  expect_equal(g$n_adjusted, 70)
  expect_error(covadj_ancova_gain(1))
})

test_that("variance factors rank the three analyses", {
  v <- covadj_variance_factors(c(0.3, 0.5, 0.7))
  expect_equal(v$change, c(1.4, 1, 0.6))
  expect_true(all(v$ancova <= pmin(v$final, v$change) + 1e-12))
})

test_that("true marginal OR is attenuated when covariate is prognostic", {
  m <- covadj_true_marginal(-0.5, 1.5, 0.7)
  expect_lt(m[["log_or"]], 0.7)
  expect_equal(covadj_true_marginal(-0.5, 0, 0.7)[["log_or"]], 0.7, tolerance = 1e-6)
})

test_that("standardize equals unadjusted with no covariates and matches arm risks", {
  d <- covadj_sim_trial(400, seed = 4)
  s <- covadj_standardize(stats::glm(y ~ trt, stats::binomial(), d), d, "trt")
  expect_equal(s$estimate[s$measure == "rd"], mean(d$y[d$trt == 1]) - mean(d$y[d$trt == 0]))
  s2 <- covadj_standardize(stats::glm(y ~ trt + x, stats::binomial(), d), d, "trt")
  expect_equal(nrow(s2), 3)
  expect_true(all(s2$lower < s2$estimate & s2$estimate < s2$upper))
})

test_that("standardized SE shrinks with a prognostic covariate and target is marginal", {
  r <- covadj_sim_compare(300, 300, seed = 5)
  tg <- covadj_true_marginal(-0.5, 1.5, 0.7)[["log_or"]]
  expect_lt(sd(r$std_est), sd(r$unadj_est))
  expect_lt(abs(mean(r$std_est) - tg), 0.05)
  expect_gt(mean(r$adj_est), tg + 0.05)
})

test_that("Bayesian g-computation: posterior SD shrinks with the covariate", {
  d <- covadj_sim_trial(300, seed = 6)
  a <- covadj_bayes_gcomp(y ~ trt, d, "trt", seed = 1)
  b <- covadj_bayes_gcomp(y ~ trt + x, d, "trt", seed = 1)
  expect_lt(sd(b$rd), sd(a$rd))
  expect_equal(nrow(a), 4000)
})
