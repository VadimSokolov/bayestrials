test_that("exact operating characteristics match Python", {
  P <- borrow_priors(30, 100)
  tabs <- lapply(P, borrow_decision_table, nt = 60, nc = 30, c = 0.975)
  expect_equal(borrow_reject_prob(tabs[[1]], .3, .3, 60, 30), 0.021872761596, tolerance = 1e-6)
  expect_equal(borrow_reject_prob(tabs[[2]], .3, .5, 60, 30), 0.787240685614, tolerance = 1e-6)
  expect_equal(borrow_reject_prob(tabs[[2]], .45, .45, 60, 30), 0.311824248374, tolerance = 1e-6)
  expect_equal(borrow_reject_prob(tabs[[3]], .3, .3, 60, 30), 0.009626053852, tolerance = 1e-6)
  expect_equal(borrow_reject_prob(tabs[[4]], .3, .5, 60, 30), 0.714211685637, tolerance = 1e-6)
})
test_that("robust weight and ESS", {
  expect_equal(borrow_robust_weight(9, 30, 30, 100), 0.9448518225579545, tolerance = 1e-6)
  wp <- borrow_robust_weight(15, 30, 30, 100)
  expect_equal(wp * 102, 70.80568936963653, tolerance = 1e-6)
  # mixture-weighted ESS: informative component has ESS nh + 2, vague Beta(1, 1) has ESS 2
  expect_equal(borrow_ess_robust(15, 30, 30, 100), wp * 102 + (1 - wp) * 2, tolerance = 1e-10)
})
test_that("basket posterior matches Python", {
  r <- borrow_basket_posterior(c(rep(12, 9), 3), rep(30, 10))
  expect_equal(r$means[1], 0.37927181, tolerance = 1e-6)
  expect_equal(r$means[10], 0.29105552, tolerance = 1e-6)
  expect_equal(r$tau_mean, 0.2679160621985298, tolerance = 1e-6)
})
test_that("MAP prior matches Python", {
  m <- borrow_map_prior(c(24, 33, 27, 45), c(80, 100, 70, 150), seed = 1)
  expect_equal(unname(quantile(m$samples, c(.025, .5, .975))), c(0.21078964, 0.32601824, 0.47205329), tolerance = 0.01)
  expect_equal(sum(m$taus * m$pw), 0.18140844553, tolerance = 1e-6)
})
test_that("robust ESS at the weight limits", {
  expect_equal(borrow_ess_robust(15, 30, 30, 100, w = 0), 2)
  expect_equal(borrow_ess_robust(15, 30, 30, 100, w = 1), 102)
})
test_that("MAP prior does not truncate heterogeneity", {
  # strongly heterogeneous history with a wide tau prior: almost all mass lies above the old cap of 2
  m <- borrow_map_prior(c(0, 100), c(100, 100), tau_scale = 5, nsamp = 2000, seed = 1)
  expect_gt(m$tau_max, 2)
  expect_lt(m$tail_mass, 1e-6)
  expect_equal(sum(m$taus * m$pw), 6.30, tolerance = 0.01)
  expect_gt(sum(m$pw[m$taus > 2]), 0.99)
  # omitted tail mass is small on an independent fine grid computed from the same model
  y <- c(0, 100); n <- c(100, 100)
  est <- log((y + .5) / (n - y + .5)); se2 <- 1 / (y + .5) + 1 / (n - y + .5)
  lp <- function(t) { w <- 1 / (se2 + t^2); mm <- sum(w * est) / sum(w)
    log(2) + dnorm(t, 0, 5, log = TRUE) + 0.5 * log(1 / sum(w)) - 0.5 * sum(log(se2 + t^2)) - 0.5 * sum(w * (est - mm)^2) }
  tg <- seq(1e-4, 400, length.out = 40000)
  d <- exp(vapply(tg, lp, 0)); d <- d / sum(d)
  expect_lt(sum(d[tg > m$tau_max]), 1e-6)
  # default homogeneous case keeps the original grid
  h <- borrow_map_prior(c(24, 33, 27, 45), c(80, 100, 70, 150), nsamp = 1000, seed = 1)
  expect_equal(h$tau_max, 2)
  expect_lt(h$tail_mass, 1e-6)
})

test_that("M3: borrow_prob_superior is exact for concentrated posteriors", {
  expect_equal(borrow_prob_superior(0, 1000, data.frame(w = 1, a = 1, b = 1001)), 0.5,
               tolerance = 1e-9)
  expect_equal(borrow_prob_superior(c(0, 500), 1000, data.frame(w = c(.3, .7), a = c(1, 1), b = c(1001, 1001)))[1], 0.5,
               tolerance = 1e-9)
  cp <- data.frame(w = c(.25, .75), a = c(10, 3), b = c(22, 9))
  yt <- c(5, 20, 41)
  ref <- vapply(yt, function(y) 0.25 * prob_superior(1 + y, 61 - y, 10, 22) +
                               0.75 * prob_superior(1 + y, 61 - y, 3, 9), numeric(1))
  expect_equal(borrow_prob_superior(yt, 60, cp), ref, tolerance = 1e-10)
  # weights are normalized
  expect_equal(borrow_prob_superior(yt, 60, transform(cp, w = 10 * w)), ref, tolerance = 1e-10)
})

test_that("M4: borrow_basket_posterior is finite for extreme data", {
  r <- borrow_basket_posterior(c(1000, 1000), c(2000, 2000))
  expect_true(all(is.finite(r$means)))
  expect_true(is.finite(r$tau_mean))
  expect_true(all(r$means > 0 & r$means < 1))
  r2 <- borrow_basket_posterior(c(0, 2000), c(2000, 2000))
  expect_true(all(is.finite(r2$means)))
  expect_lt(r2$means[1], r2$means[2])
  # unchanged for moderate data: symmetric baskets give equal means
  r3 <- borrow_basket_posterior(c(12, 12), c(30, 30), nq = 21)
  expect_equal(r3$means[1], r3$means[2], tolerance = 1e-12)
})
