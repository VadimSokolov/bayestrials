test_that("exact helpers match Python", {
  tb <- platform_shared_control_table(2:3, 80)
  expect_equal(tb$platform, c(240, 320))
  expect_equal(tb$saving, c(0.25, 1/3), tolerance = 1e-12)
  expect_equal(platform_contrast_variance(4, 480, 2), 0.01875, tolerance = 1e-12)
  expect_equal(platform_arm_correlation(1), 0.5)
})
test_that("false positive counts", {
  r <- platform_false_positive_counts(3, 1, seed = 1)
  expect_equal(sum(r), 1)
  se <- sqrt(c(0.9374075, 0.0516125, 0.00932, 0.00166) * (1 - c(0.9374075, 0.0516125, 0.00932, 0.00166)) / 4e5)
  expect_true(all(abs(r - c(0.9374075, 0.0516125, 0.00932, 0.00166)) < 4 * se * sqrt(2)))
})
test_that("drift simulation matches Python", {
  r <- platform_drift_simulation(nsim = 4000, seed = 7)
  ref_rej <- c(0.02625, 0.08925, 0.028)
  se <- sqrt(ref_rej * (1 - ref_rej) / 4000)
  expect_true(all(abs(r$rejection - ref_rej) < 4 * se * sqrt(2)))
  expect_equal(r$mean_estimate[2], 0.02472, tolerance = 0.2)
})
test_that("predictive graduation matches Python", {
  r <- platform_predictive_graduation(20, 40, 12, 40, seed = 3)
  expect_lt(abs(r[["pp"]] - 0.761525), 4 * sqrt(0.76 * 0.24 / 4e4) * sqrt(2))
  expect_lt(abs(r[["post"]] - 0.96515), 4 * sqrt(0.97 * 0.03 / 4e4) * sqrt(2))
})
test_that("time machine summaries are sensible", {
  r <- platform_time_machine(c(0, .02, .04, .06), seed = 11)
  expect_equal(r$method[1:3], c("current only", "pool all", "learned smoothing (REML tau)"))
  expect_equal(r$bias[2], -0.0298, tolerance = 0.05)  # pooling bias is exactly mean(shift) - last shift
  expect_equal(r$rmse[1], 0.04775, tolerance = 0.05)  # binomial sd of the last period mean
  expect_true(all(is.finite(r$rmse)))
  # learned smoothing is competitive with the better of the two extreme estimators
  expect_lt(r$rmse[3], max(r$rmse[1:2]))
  # under a stable control rate smoothing beats the current-only estimate; under a jump it beats pooling
  s0 <- platform_time_machine(rep(0, 4), nsim = 2000, seed = 2)
  expect_lt(s0$rmse[3], s0$rmse[1])
  s1 <- platform_time_machine(c(0, 0, .15, .15), nsim = 2000, seed = 2)
  expect_lt(s1$rmse[3], s1$rmse[2])
})
test_that("time machine fit uses only observed counts, not the true rates", {
  fit <- bayestrials:::.tm_fit
  expect_false(any(c("truth", "p0", "shift") %in% names(formals(fit))))
  y <- rbind(c(30, 28, 35, 33), c(10, 12, 40, 41))
  a <- fit(y, 100); b <- fit(y, 100)
  expect_identical(a, b)
  # identical counts, however they were generated, give identical estimates
  set.seed(1); y1 <- matrix(rbinom(8, 100, 0.3), 2)
  set.seed(1); y2 <- matrix(rbinom(8, 100, rep(0.3, 8)), 2)
  expect_identical(fit(y1, 100), fit(y2, 100))
  # smoothing is learned: stable counts give a small tau, a jump gives a larger tau and a weight near current
  stable <- fit(rbind(c(30, 30, 30, 30)), 100)
  jump <- fit(rbind(c(10, 10, 45, 45)), 100)
  expect_lt(stable$tau_hat, jump$tau_hat)
  expect_gt(jump$est[[1]], 0.35)
  expect_lt(abs(stable$est[[1]] - 0.30), 0.005)
})
