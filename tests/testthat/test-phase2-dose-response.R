test_that("simon_optimal reproduces the published table and exact error rates", {
  s <- simon_optimal(0.2, 0.4, 0.05, 0.80)
  expect_equal(c(s$r1, s$n1, s$r, s$n), c(3, 13, 12, 43))
  expect_equal(s$EN, 20.58, tolerance = 1e-3)
  expect_equal(simon_oc(s$n1, s$r1, s$n, s$r, 0.2)[["P_reject"]], 0.0496, tolerance = 1e-2)
  expect_equal(simon_oc(s$n1, s$r1, s$n, s$r, 0.4)[["P_reject"]], 0.800, tolerance = 1e-2)
})

test_that("mcp_contrasts critical value and mcp_power null behavior", {
  skip_if_not_installed("mvtnorm")
  d <- c(0, 10, 25, 50, 100)
  sh <- cbind(linear = d / 100, emax = d / (d + 15), exponential = (exp(d / 50) - 1) / (exp(2) - 1),
              umbrella = 1 - (d / 50 - 1)^2)
  m <- mcp_contrasts(sh)
  expect_equal(colSums(m$contrasts^2), c(linear = 1, emax = 1, exponential = 1, umbrella = 1))
  expect_gt(m$crit, qnorm(0.975)); expect_lt(m$crit, qnorm(1 - 0.025 / 4))
  p <- mcp_power(rep(0, 5), m$contrasts, m$crit, n = 40, B = 20000, seed = 1)
  expect_equal(unname(p[["mcp"]]), 0.025, tolerance = 0.2)
})
