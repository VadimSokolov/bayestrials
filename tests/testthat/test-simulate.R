test_that("post_table matches Python", {
  expect_equal(oc_post_table(10, 0)[4, 8], 0.956945531248932, tolerance = 1e-6)
  expect_equal(oc_post_table(10, 0.05)[4, 8], 0.9293309606529938, tolerance = 1e-6)
})

test_that("predictive_table matches Python", {
  expect_equal(oc_predictive_table(5, 10, 0.9)[3, 4], 0.35494555949101403, tolerance = 1e-6)
  expect_equal(oc_predictive_table(10, 10, 0.9)[3, 4], 0, tolerance = 1e-6)
})

test_that("calibration and operating characteristics agree with Python within MC error", {
  looks <- c(50, 100, 150, 200)
  nsim <- 20000
  s0 <- oc_simulate(0.3, 0.3, looks, nsim = nsim, seed = 1)
  cc <- oc_calibrate(s0)
  expect_lt(abs(cc - 0.9906811), 0.003)
  expect_lte(oc_summarize(s0, cc)$reject, 0.025)
  s1 <- oc_simulate(0.3, 0.5, looks, nsim = nsim, seed = 2, fut = 0.1)
  o <- oc_summarize(s1, 0.9906811)
  expect_lt(abs(o$reject - 0.969), 4 * sqrt(2) * o$se_reject)
  expect_lt(abs(o$EN - 198.48), 3)
  expect_lt(abs(o$bias_succ - 0.0324), 0.01)
})

test_that("normal method runs and futility masks later looks", {
  s <- oc_simulate(0.3, 0.3, c(50, 100), nsim = 500, seed = 3, method = "normal", fut = 0.2)
  expect_true(all(is.na(s$pp[!s$alive])))
  expect_equal(dim(oc_predictive_matrix(oc_simulate(0.3, 0.3, c(5, 10), nsim = 20, seed = 1), 10, 0.9)),
               c(20, 2))
})
