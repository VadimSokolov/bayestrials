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
  expect_equal(borrow_ess_robust(15, 30, 30, 100), 70.80568936963653, tolerance = 1e-6)
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
