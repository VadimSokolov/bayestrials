test_that("rofecoxib keeps supplied data", {
  r <- safety_rofecoxib(rr = .5, lo = .4, hi = .6, period = "x", s0 = Inf)
  expect_equal(r$post_rr, 0.5)
})
test_that("decide_ss_min_n does not assume monotonicity", {
  expect_equal(decide_ss_min_n(.99, analysis_mean = .3, analysis_sd = .1), 1)
})
test_that("reanalysis_hier respects tau_max", {
  h <- reanalysis_hier(c(.00044, .0002, -.00017), c(.00012, .00007, .00013), tau_max = .0001)
  expect_equal(h$tau_mean, 6.5334e-05, tolerance = 1e-3)
  expect_lte(h$tau_mean, 1e-4)
  expect_equal(h$strata$p_benefit[3], 0.8664, tolerance = 1e-3)
})
test_that("surv_gamma_hazard median mean is Inf without a finite mean", {
  expect_identical(unname(surv_gamma_hazard(0, 100)$median["mean"]), Inf)
})
test_that("GPS fit likelihood does not underflow", {
  set.seed(1)
  f <- suppressWarnings(safety_gps_fit(c(9e5, 1e6, 1.1e6), rep(1, 3)))
  expect_true(all(is.finite(f$par)))
  # with underflow the objective was flat and the optimizer stayed near the start (loglik about -2072 vs -285973 true)
  mu <- c(f$par[["a1"]] / f$par[["b1"]], f$par[["a2"]] / f$par[["b2"]])
  expect_true(any(abs(mu - 1e6) < 2e5))
  expect_gt(f$loglik, -100)
})
test_that("monitor_calibrate_bayes brackets small roots", {
  m <- suppressWarnings(monitor_calibrate_bayes(100, prior_sd = .1))
  expect_equal(m$threshold, 0.809626986, tolerance = 1e-6)
  expect_equal(m$bounds, 1.959963985, tolerance = 1e-6)
})
test_that("safety_gps quantile above 1000", {
  expect_equal(safety_gps(200, .01)$eb05, 1613.71684, tolerance = 1e-6)
})
test_that("approximate bar honors c", {
  t <- decide_bayes_type1_table(n = 100, c = 10)
  expect_equal(t$bar_approx, 3.03649, tolerance = 1e-5)
  expect_equal(t$alpha_approx, 0.00239348, tolerance = 1e-4)
})
test_that("stable Bayes factors", {
  expect_equal(decide_bf_q(40, 1e100), 2.72637e297, tolerance = 1e-5)
  expect_equal(safety_bf_monitor(est = 40, se = 1, tau = .01)$bf, 1.08322424, tolerance = 1e-6)
})
test_that("jeffreys exact recycles and centers", {
  expect_equal(decide_jeffreys_k_exact(2, c(10, 100)), c(.75285095, 1.77821860), tolerance = 1e-6)
  expect_equal(decide_jeffreys_k_exact(35, 10), 4.82292e-264, tolerance = 1e-4)
})
test_that("surv_pp_cp final look", {
  r <- surv_pp_cp(600, .8)
  expect_equal(c(r$pp, r$cp_design, r$cp_trend), c(1, 1, 1))
  expect_error(surv_pp_cp(700, .8))
})
test_that("single simulated look works", {
  s <- safety_bf_monitor(looks = 1, nsim = 10)
  expect_length(s$false_alarm, 3)
})
test_that("decide_bayes_bar edge cases", {
  expect_equal(decide_bayes_bar(1, c = .5), 0)
  expect_equal(decide_bayes_bar(1e-20), 1, tolerance = 1e-6)
})
test_that("decide_crossing_q adaptive", {
  expect_equal(decide_crossing_q(1.1), 0.483155037, tolerance = 1e-6)
  expect_equal(decide_crossing_q(1.96), 41.6, tolerance = 0.01)
  expect_error(decide_crossing_q(0.9))
})
test_that("surv_rmst near zero hazard", {
  expect_equal(surv_rmst(1, 1e-17), 1, tolerance = 1e-8)
})
test_that("decide_inb handles negative conversion", {
  r <- decide_inb(-.22, .058, -.25, 1e5, 4000)
  expect_equal(r$sd, 1450, tolerance = 1e-8)
  expect_equal(r$prob_positive, 0.849545, tolerance = 1e-5)
})
test_that("BLRM grid is accurate at the reference dose", {
  t <- dose_blrm_table()
  got <- unlist(t[4, c("p_under", "p_target", "p_over")], use.names = FALSE)
  expect_lt(max(abs(got - c(0.359851, 0.450561, 0.189589))), 5e-4)
})
