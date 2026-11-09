test_that("core boundaries match published values", {
  expect_equal(gs_boundary(4), c(4.049, 2.863, 2.337, 2.024), tolerance = 1e-3)
  expect_equal(gs_boundary(4, type = "pocock"), rep(2.361, 4), tolerance = 1e-3)
  expect_equal(gs_spending_boundary((1:4) / 4), c(4.333, 2.963, 2.359, 2.014), tolerance = 1e-3)
})

test_that("two-sided recursion matches Python", {
  t <- (1:4) / 4
  b <- gs_boundary(4)
  expect_equal(gs_power(t, b), 0.025, tolerance = 1e-6)
  expect_equal(gs_power_drift(t, b), 3.2772396, tolerance = 1e-5)
  lo <- c(-9, 0, 0.5, 1.0)
  expect_equal(gs_power_drift(t, b, lo), 3.2898853, tolerance = 1e-5)
  expect_equal(gs_expected_information(t, b, lo, 0), 0.6920757, tolerance = 1e-5)
  expect_equal(gs_expected_information(t, b, lo, 3), 0.7782860, tolerance = 1e-5)
  expect_equal(gs_expected_information(t, b, NULL, 3), 0.7893430, tolerance = 1e-5)
  r <- gs_crossing_two_sided(t, b, lo, 3)
  expect_equal(r$up, c(0.00540795, 0.22401125, 0.37834727, 0.23227054), tolerance = 1e-6)
  expect_equal(r$down[2:4], c(1.69474267e-02, 1.03675129e-02, 1.09818068e-02), tolerance = 1e-6)
  r0 <- gs_crossing_two_sided(t, rep(1.96, 4))
  expect_equal(r0$up, c(0.0249979, 0.01655768, 0.01206932, 0.00945994), tolerance = 1e-6)
  expect_equal(gs_crossing_prob(rep(1.96, 4), t)$cross, r0$up, tolerance = 1e-5)
})

test_that("naive inflation and tolls", {
  expect_lt(abs(gs_naive_inflation(4) - 0.0638), 2e-3)
  expect_equal(gs_naive_inflation(1), 0.025, tolerance = 1e-3)
  expect_lt(abs(gs_naive_inflation(2) - 0.042), 2e-3)
  tl <- gs_tolls(c(25, 4000), tau = 0.25)
  expect_equal(tl$surcharge[1], 7.378, tolerance = 1e-4)
  expect_equal(tl$ockham, log(1 + c(25, 4000) * 0.0625))
  expect_equal(tl$anytime_bar^2 - tl$bayes_bar^2, rep(tl$surcharge[1], 2))
})

test_that("gs_tolls offers exact finite-information bars next to the asymptotic ones", {
  tl <- gs_tolls(c(1, 25, 4000), tau = 0.25)
  q <- c(1, 25, 4000) * 0.0625
  expect_equal(tl$bayes_bar_exact, sqrt((1 + 1 / q) * log(1 + q)), tolerance = 1e-12)
  expect_equal(tl$anytime_bar_exact, sqrt((1 + 1 / q) * (log(1 + q) + 2 * log(1 / 0.025))),
               tolerance = 1e-12)
  # the exact bars are exactly where the mixture Bayes factor hits 1 and 1/alpha
  expect_equal(bf10_normal(tl$bayes_bar_exact, q), rep(1, 3), tolerance = 1e-10)
  expect_equal(bf10_normal(tl$anytime_bar_exact, q), rep(40, 3), tolerance = 1e-10)
  # finite-information factor matters at small q (review: 0.246 vs about 1.015)
  expect_equal(tl$bayes_bar[1], 0.2462, tolerance = 1e-3)
  expect_equal(tl$bayes_bar_exact[1], 1.0152, tolerance = 1e-3)
  expect_equal(tl$ratio_exact, tl$anytime_bar_exact^2 / tl$bayes_bar_exact^2)
  # the two forms agree as q grows
  big <- gs_tolls(1e8, tau = 0.25)
  expect_equal(big$bayes_bar_exact, big$bayes_bar, tolerance = 1e-6)
})

mc_cross <- function(t, b, nrep = 2e5, drift = 0, lower = NULL, seed = 1) {
  set.seed(seed)
  K <- length(t); d <- diff(c(0, t))
  S <- matrix(rnorm(nrep * K, rep(drift * d, each = nrep), rep(sqrt(d), each = nrep)), nrep, K)
  for (j in seq_len(K)[-1]) S[, j] <- S[, j] + S[, j - 1]
  Z <- sweep(S, 2, sqrt(t), "/")
  alive <- rep(TRUE, nrep); up <- numeric(K); dn <- numeric(K)
  for (k in seq_len(K)) {
    u <- alive & Z[, k] > b[k]; up[k] <- mean(u)
    if (!is.null(lower)) { l <- alive & Z[, k] < lower[k]; dn[k] <- mean(l); alive <- alive & !l }
    alive <- alive & !u
  }
  list(up = up, down = dn)
}

test_that("M1: nearly coincident and irregular looks give valid probabilities", {
  two <- gs_crossing_prob(rep(1.96, 2), c(0.5, 1))$total
  expect_equal(two, 0.04156, tolerance = 1e-3)
  for (d in c(1e-2, 1e-4, 1e-6)) {
    r <- gs_crossing_prob(rep(1.96, 3), c(0.5, 0.5 + d, 1))
    expect_lte(r$total, 1)
    expect_true(all(r$cross >= 0))
    # continuity: an almost coincident look adds at most O(sqrt(d)) crossing mass
    expect_lt(abs(r$total - two), 0.3 * sqrt(d) + 1e-5)
  }
  # failing input from the review (previously about 1.85)
  expect_lt(gs_crossing_prob(rep(1.96, 3), c(0.5, 0.50000001, 1))$total, 0.0416 + 1e-3)
  # irregular spacing against Monte Carlo
  t <- c(0.05, 0.0501, 0.4, 0.9, 1); b <- c(2.6, 2.6, 2.4, 2.1, 2.0)
  r <- gs_crossing_two_sided(t, b, drift = 1)
  m <- mc_cross(t, b, drift = 1)
  expect_true(all(abs(r$up - m$up) < 4 * sqrt(m$up * (1 - m$up) / 2e5) + 1e-4))
  # quadrature convergence in the node count
  r2 <- gs_crossing_two_sided(t, b, drift = 1, nnode = 600)
  expect_equal(r$up, r2$up, tolerance = 1e-7)
  # probability mass conservation with a futility boundary
  tt <- c(0.1, 0.3, 0.31, 1)
  g <- gs_crossing_two_sided(tt, c(3, 2.5, 2.4, 2), lower = c(-1, 0, 0.2, 1), drift = 2)
  expect_lte(sum(g$up) + sum(g$down), 1)
})

test_that("M2: supplied futility boundaries are applied exactly", {
  expect_equal(gs_crossing_two_sided(1, 2, lower = -10, drift = -10)$down, 0.5)
  # two looks, boundary below the old -9 sqrt(t) clamp: closed form at look 1
  t <- c(0.5, 1)
  r <- gs_crossing_two_sided(t, c(10, 10), lower = c(-10, -10), drift = -10)
  expect_equal(r$down[1], pnorm((-10 * sqrt(0.5) + 10 * 0.5) / sqrt(0.5)), tolerance = 1e-12)
  expect_equal(r$up, c(0, 0), tolerance = 1e-12)
  expect_error(gs_crossing_two_sided(c(0.5, 1), c(1, 1), lower = c(2, 0)))
})

test_that("input validation for looks and boundaries", {
  expect_error(gs_crossing_prob(c(2, 2), c(0.5, 0.4)))
  expect_error(gs_crossing_prob(c(2, 2), c(0, 1)))
  expect_error(gs_crossing_prob(c(2, 2, 2), c(0.5, 1)))
  expect_error(gs_boundary(4, alpha = 0))
  expect_error(gs_spending_boundary(c(0.5, 0.5)))
})

test_that("MINOR 2: root brackets expand adaptively", {
  expect_equal(gs_boundary(1, alpha = 0.2), qnorm(0.8), tolerance = 1e-6)
  b <- gs_boundary(4, alpha = 0.4)
  expect_equal(gs_crossing_prob(b, (1:4) / 4)$total, 0.4, tolerance = 1e-8)
  b2 <- gs_boundary(3, alpha = 1e-9, type = "pocock")
  expect_equal(gs_crossing_prob(b2, (1:3) / 3)$total, 1e-9, tolerance = 1e-5)
  t <- (1:4) / 4
  bb <- gs_boundary(4)
  d <- gs_power_drift(t, bb, target = 0.03)
  expect_equal(gs_power(t, bb, drift = d), 0.03, tolerance = 1e-8)
  d2 <- gs_power_drift(t, bb, target = 0.999)
  expect_equal(gs_power(t, bb, drift = d2), 0.999, tolerance = 1e-8)
})

test_that("MINOR 1: Lan-DeMets spending is stable for negligible early spending", {
  b <- gs_spending_boundary(c(0.001, 0.5, 1))
  expect_true(all(is.finite(b)))
  expect_gt(b[1], 10)
  expect_equal(b[2:3], gs_spending_boundary(c(0.5, 1)), tolerance = 1e-4)
  expect_equal(gs_crossing_prob(b, c(0.001, 0.5, 1))$total, 0.025, tolerance = 1e-6)
  # boundaries reproduce the spending increments
  t <- c(0.2, 0.45, 0.8, 1)
  bs <- gs_spending_boundary(t)
  cum <- 2 * pnorm(qnorm(0.0125, lower.tail = FALSE) / sqrt(t), lower.tail = FALSE)
  expect_equal(cumsum(gs_crossing_prob(bs, t)$cross), cum, tolerance = 1e-6)
  bp <- gs_spending_boundary(t, spending = "ld_pocock")
  expect_equal(sum(gs_crossing_prob(bp, t)$cross), 0.025, tolerance = 1e-6)
  # nearly coincident looks do not break the spending recursion
  bc <- gs_spending_boundary(c(0.5, 0.5 + 1e-6, 1))
  expect_true(all(is.finite(bc)))
  expect_equal(gs_crossing_prob(bc, c(0.5, 0.5 + 1e-6, 1))$total, 0.025, tolerance = 1e-6)
})
