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

test_that("M3: oc_post_table is exact on concentrated and U-shaped posteriors", {
  # exchangeable states give exactly 1/2 (the old fixed grid gave 0.377 at n = 1000, Beta(.5, 1000.5))
  expect_equal(prob_superior(.5, 1000.5, .5, 1000.5), 0.5, tolerance = 1e-9)
  tb <- oc_post_table(120, prior = c(.5, .5))
  expect_equal(tb[1, 1], 0.5, tolerance = 1e-9)
  expect_equal(tb[121, 121], 0.5, tolerance = 1e-9)
  expect_equal(tb + t(tb), matrix(1, 121, 121), tolerance = 1e-9)
  tj <- oc_post_table(6, prior = c(.01, .01))
  expect_equal(diag(tj), rep(0.5, 7), tolerance = 1e-9)
  # agrees with the cell-by-cell exact integral
  tt <- oc_post_table(12, a = 0.1, prior = c(.5, .5))
  expect_equal(tt[3, 9], prob_superior(.5 + 8, .5 + 4, .5 + 2, .5 + 10, margin = 0.1), tolerance = 1e-9)
  # complementary probabilities across margins: Pr(p1 - p0 > a) + Pr(p0 - p1 > -a) = 1
  ta <- oc_post_table(8, a = 0.2); tn <- oc_post_table(8, a = -0.2)
  expect_equal(ta + t(tn), matrix(1, 9, 9), tolerance = 1e-8)
})

test_that("MINOR 9: posterior-table cache keys are exact", {
  t0 <- oc_post_table(10, a = 0)
  t1 <- oc_post_table(10, a = 4e-7)
  expect_false(identical(t0, t1))
  expect_gt(max(abs(t1 - t0)), 0)
  tp <- oc_post_table(10, prior = c(1, 1))
  tq <- oc_post_table(10, prior = c(1, 1 + 1e-9))
  expect_false(identical(tp, tq))
  expect_identical(oc_post_table(10, a = 4e-7), t1)
})

test_that("MINOR 12: calibration is an empirical statement about the supplied null", {
  looks <- c(40, 80)
  n1 <- oc_simulate(0.3, 0.3, looks, nsim = 4000, seed = 1)
  cc <- oc_calibrate(n1)
  expect_lte(oc_summarize(n1, cc)$reject, 0.025)
  fresh <- oc_simulate(0.3, 0.3, looks, nsim = 4000, seed = 2)
  expect_lt(abs(oc_summarize(fresh, cc)$reject - 0.025), 0.012)
  rd <- tools::Rd_db("bayestrials")
  if (length(rd)) {
    txt <- paste(capture.output(print(rd[["oc_calibrate.Rd"]])), collapse = " ")
    expect_match(txt, "empirical", ignore.case = TRUE)
    expect_match(txt, "nuisance")
    expect_match(txt, "Monte Carlo")
    expect_match(txt, "fresh")
  }
})
