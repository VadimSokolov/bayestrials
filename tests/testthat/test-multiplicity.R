if (!exists("ulcer")) load(testthat::test_path("..", "..", "data", "ulcer.rda"), envir = globalenv())

test_that("step procedures match Python", {
  p <- c(0.001, 0.008, 0.012, 0.4)
  expect_equal(mult_bonferroni(p), c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(mult_holm(p), c(TRUE, TRUE, TRUE, FALSE))
  expect_equal(mult_hochberg(c(0.011, 0.012, 0.5)), c(TRUE, TRUE, FALSE))
  expect_equal(mult_holm(c(0.011, 0.012, 0.5)), c(FALSE, FALSE, FALSE))
})

test_that("graphical procedure passes alpha along the graph", {
  G <- matrix(c(0, 1, 1, 0), 2)
  expect_equal(mult_graph(c(0.04, 0.001), c(0.5, 0.5), G, alpha = 0.05), c(TRUE, TRUE))
  expect_equal(mult_graph(c(0.04, 0.03), c(0.5, 0.5), G, alpha = 0.05), c(FALSE, FALSE))
})

test_that("ulcer numbers reproduce chapter 9", {
  l <- shrink_logodds()
  d <- l$y[1] - l$y[4]; se <- sqrt(l$v[1] + l$v[4])
  p2 <- 2 * (1 - pnorm(d / se))
  expect_equal(round(min(1, 6 * p2), 4), 0.0855)
  fit <- shrink_eb_fit(l$y, l$v)
  expect_equal(fit$tau, 0.26441641439, tolerance = 1e-4)
  expect_equal(fit$mu, -1.8291455616, tolerance = 1e-5)
  expect_equal(fit$shrunk[1], -1.56952578, tolerance = 1e-5)
  dsh <- fit$shrunk[1] - fit$shrunk[4]; ssh <- sqrt(fit$sd[1]^2 + fit$sd[4]^2)
  expect_equal(round(pnorm(dsh / ssh), 3), 0.952)
})

test_that("simulation of procedures is sensible", {
  r <- mult_simulate(4, c(0, 0, 0, 0), nsim = 4000, seed = 1)
  expect_gt(r$fwer[1], 0.05)
  expect_lt(r$fwer[2], 0.04)
  expect_true(all(r$fwer[3:4] <= r$fwer[1]))
})
