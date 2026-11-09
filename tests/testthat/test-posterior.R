test_that("hpd_beta matches Python", {
  h <- hpd_beta(1.25 + 3, 5.75 + 10)
  expect_equal(unname(h), c(0.05363968784829362, 0.3880881025663918), tolerance = 1e-6)
  h2 <- hpd_beta(8, 4, 0.9)
  expect_equal(unname(h2), c(0.45844636507428854, 0.8825246721408506), tolerance = 1e-6)
})

test_that("sim_sequential_looks agrees with Python within MC error", {
  r <- sim_sequential_looks(4000, 50, 5, seed = 11)
  se_bin <- function(p, n) sqrt(p * (1 - p) / n)
  expect_lt(abs(r$final_only - 0.02725), 4 * se_bin(0.027, 4000))
  expect_lt(abs(r$any_look - 0.0705), 4 * se_bin(0.07, 4000))
  expect_lt(abs(r$n_stopped - 1732), 4 * sqrt(4000 * 0.43 * 0.57))
  expect_lt(abs(r$mean_post_at_stop - 0.9965), 0.002)
  expect_lt(abs(r$frac_truly_better - 0.99596), 0.004)
})

test_that("dbetabinom matches choose/beta formula", {
  x <- 0:50
  expect_equal(dbetabinom(x, 50, 46, 56), choose(50, x) * beta(x + 46, 50 - x + 56) / beta(46, 56), tolerance = 1e-10)
  expect_equal(sum(dbetabinom(x, 50, 46, 56)), 1)
  expect_equal(dbetabinom(3, 10, 2, 3, log = TRUE), log(dbetabinom(3, 10, 2, 3)))
  expect_equal(dbetabinom(c(-1, 11), 10, 2, 3), c(0, 0))
})

test_that("grid_update matches conjugate Beta posterior in the limit", {
  theta <- seq(0.0005, 0.9995, by = 0.001)
  g <- grid_update(theta, dbeta(theta, 2, 2), dbinom(7, 10, theta))
  expect_equal(sum(g$posterior), 1)
  expect_equal(sum(g$theta * g$posterior), 9 / 14, tolerance = 1e-4)
  expect_equal(sum(g$prior), 1)
  expect_error(grid_update(c(0, 1), c(1, 0), c(0, 1)))
})

test_that("prob_superior exact: exchangeability, complements, closed forms", {
  expect_equal(prob_superior(5, 7, 5, 7), 0.5, tolerance = 1e-9)
  expect_equal(prob_superior(0.5, 0.5, 0.5, 0.5), 0.5, tolerance = 1e-9)   # Jeffreys
  expect_equal(prob_superior(1, 1, 1, 1), 0.5, tolerance = 1e-9)
  # review failures: U-shaped and concentrated identical Betas
  expect_equal(prob_superior(0.01, 0.01, 0.01, 0.01), 0.5, tolerance = 1e-6)
  expect_equal(prob_superior(1e6, 1, 1e6, 1), 0.5, tolerance = 1e-6)
  expect_equal(prob_superior(1, 1e6, 1, 1e6), 0.5, tolerance = 1e-6)
  expect_equal(prob_superior(1e5, 1e5, 1e5, 1e5), 0.5, tolerance = 1e-8)
  # complementary probabilities
  a <- c(2.5, 0.3, 40); b <- c(7, 0.9, 3); c0 <- c(3, 1.2, 35); d0 <- c(4, 0.2, 5)
  expect_equal(prob_superior(a, b, c0, d0) + prob_superior(c0, d0, a, b), rep(1, 3),
               tolerance = 1e-8)
  # uniform control: Pr(p1 > p0 + m) = E[(p1 - m)+]
  expect_equal(prob_superior(3, 2, 1, 1), 0.6, tolerance = 1e-9)
  expect_equal(prob_superior(1, 1, 1, 1, margin = 0.2), 0.8^2 / 2, tolerance = 1e-9)
  expect_equal(prob_superior(1, 1, 1, 1, margin = -0.2), 1 - 0.8^2 / 2, tolerance = 1e-9)
  # integer shapes: Cook's finite-sum formula
  cook <- function(a, b, c, d) sum(sapply(0:(a - 1), function(j)
    beta(c + j, b + d) / ((b + j) * beta(1 + j, b) * beta(c, d))))
  expect_equal(prob_superior(7, 5, 4, 9), cook(7, 5, 4, 9), tolerance = 1e-9)
  expect_equal(prob_superior(3, 11, 6, 6), cook(3, 11, 6, 6), tolerance = 1e-9)
  # endpoint counts: 0/1000 versus 0/1000 with a Jeffreys prior, and clear separation
  expect_equal(prob_superior(.5, 1000.5, .5, 1000.5), 0.5, tolerance = 1e-8)
  expect_equal(prob_superior(1e6, 1, 1, 1e6), 1, tolerance = 1e-9)
  expect_equal(prob_superior(1, 1e6, 1e6, 1), 0, tolerance = 1e-9)
  # concentrated versus diffuse: arm 1 is essentially the point 0.5
  expect_equal(prob_superior(1e7, 1e7, 2, 2), 0.5, tolerance = 1e-5)
  # margins that clip
  expect_equal(prob_superior(3, 2, 1, 1, margin = 1), 0)
  expect_equal(prob_superior(3, 2, 1, 1, margin = -1), 1)
  expect_equal(prob_superior(3, 2, 1, 1, margin = 0.9),
               integrate(function(p) (p - 0.9) * dbeta(p, 3, 2), 0.9, 1)$value, tolerance = 1e-9)
  expect_error(prob_superior(0, 1, 1, 1))
  expect_error(prob_superior(1, NA, 1, 1))
  expect_error(prob_superior(1, 1, 1, 1, margin = Inf))
})

test_that("prob_superior exact agrees with Monte Carlo and the normal method", {
  set.seed(5)
  n <- 4e5
  cases <- list(c(0.2, 0.3, 2, 0.5, 0), c(12, 30, 8, 33, 0), c(4, 6, 3, 7, 0.1), c(0.5, 0.5, 20, 20, -0.4))
  for (cs in cases) {
    ex <- prob_superior(cs[1], cs[2], cs[3], cs[4], cs[5])
    mc <- mean(rbeta(n, cs[1], cs[2]) > rbeta(n, cs[3], cs[4]) + cs[5])
    expect_lt(abs(ex - mc), 4 * sqrt(0.25 / n))
  }
  mcm <- prob_superior(12, 30, 8, 33, method = "mc", ndraw = 2e5)
  expect_lt(abs(mcm - prob_superior(12, 30, 8, 33)), 4 * sqrt(0.25 / 2e5))
  expect_equal(prob_superior(120, 300, 80, 330, method = "normal"),
               prob_superior(120, 300, 80, 330), tolerance = 5e-3)
})

test_that("beta_post validates counts and prior shapes (M14)", {
  expect_equal(beta_post(3, 10), list(a = 4, b = 8))
  expect_error(beta_post(2, 1))
  expect_error(beta_post(-1, 5))
  expect_error(beta_post(1.5, 5))
  expect_error(beta_post(1, 5, a = 0))
  expect_error(beta_post(1, 5, b = -1))
})

test_that("sim_sequential_looks keeps matrix shape with a single look or trial (MINOR 8)", {
  r1 <- sim_sequential_looks(200, 20, 1, draws = 30, seed = 1)
  expect_true(is.finite(r1$final_only))
  expect_equal(r1$final_only, r1$any_look)
  r2 <- sim_sequential_looks(1, 20, 3, draws = 30, seed = 2)
  expect_true(r2$n_stopped %in% 0:1)
  r3 <- sim_sequential_looks(1, 20, 1, draws = 30, seed = 3)
  expect_true(r3$n_stopped %in% 0:1)
  expect_error(sim_sequential_looks(10, 20, 0))
})

test_that("hpd_beta warns outside the unimodal range", {
  expect_warning(hpd_beta(0.5, 2))
  expect_error(hpd_beta(-1, 2))
})

test_that("N1: prob_superior resolves reflected-quantile endpoint mass", {
  expect_equal(prob_superior(.015, 1, .03, 1), 1 / 3, tolerance = 1e-8)
  expect_equal(prob_superior(.03, 1, .015, 1), 2 / 3, tolerance = 1e-8)
  expect_equal(prob_superior(.05, 1, .1, 1), 1 / 3, tolerance = 1e-8)
  expect_equal(prob_superior(1, .015, 1, .03), 2 / 3, tolerance = 1e-8)
  expect_equal(prob_superior(1, 1, 1, 1, margin = 0.2), 0.8^2 / 2, tolerance = 1e-9)
})
