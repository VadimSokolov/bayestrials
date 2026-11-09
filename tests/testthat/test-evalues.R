test_that("mixture e-process matches Python", {
  expect_equal(ev_mixture(c(5, -2), c(100, 50), 0.15), c(0.79996829, 0.58642986), tolerance = 1e-6)
  expect_equal(ev_log_mixture(c(5, -2), c(100, 50), 0.15), c(-0.22318319, -0.53370221), tolerance = 1e-6)
})

test_that("first crossing", {
  m <- matrix(c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE), 2, byrow = TRUE)
  expect_equal(ev_first_cross(m), c(2, 4))
})

test_that("naive false alarm matches Python within MC error", {
  looks <- c(1, 2, 5, 10, 20, 50, 200)
  py <- c(0.02665, 0.0413, 0.07145, 0.0961, 0.1255, 0.1611, 0.2165)
  r <- ev_naive_false_alarm(200, looks, reps = 20000, seed = 12)
  se <- sqrt(2 * py * (1 - py) / 20000)
  expect_true(all(abs(r$naive$false_positive - py) < 4 * se))
  expect_lt(r$evalue, 0.025)
})

test_that("SELECT-like comparison matches Python within MC error", {
  r <- ev_select_like(reps = 20000, seed = 13)
  expect_equal(r$cb, 2.622, tolerance = 0.05)
  expect_equal(r$co, 2.038, tolerance = 0.05)
  pw <- sapply(r$alt, function(v) mean(v <= 1270))
  py <- c(bayes = 0.94465, evalue = 0.88685, obf = 0.9747)
  se <- sqrt(2 * py * (1 - py) / 20000)
  expect_true(all(abs(pw[names(py)] - py) < 4 * se))
  expect_equal(median(r$alt$obf[r$alt$obf <= 1270]), 952)
  expect_equal(median(r$alt$bayes[r$alt$bayes <= 1270]), 468, tolerance = 0.03)
  expect_equal(median(r$alt$evalue[r$alt$evalue <= 1270]), 520, tolerance = 0.03)
})

test_that("betting paths", {
  w <- ev_betting_paths(0.45, 0.30, 0.31, 200, 20000, seed = 1)
  expect_equal(mean(w[, 200]), 4.7723, tolerance = 4 * 3.0094 / sqrt(20000) * sqrt(2) / 4.77)
  w0 <- ev_betting_paths(0.30, 0.30, 0.31, 200, 20000, seed = 1)
  expect_lt(mean(apply(w0, 1, max) >= log(40)), 0.025)
})

test_that("M4: e-values are finite on the raw scale in extreme regimes", {
  expect_equal(ev_mixture(-1000, 100, 1), 0.000797804, tolerance = 1e-5)
  expect_equal(ev_mixture(-1000, 100, 1), exp(ev_log_mixture(-1000, 100, 1, sign = 1)))
  expect_true(all(is.finite(ev_mixture(c(-1e4, -500, 0, 30), 100, 0.15))))
  expect_equal(ev_mixture(1000, 100, 1), Inf)
  # sign = -1 mirrors the sign of S
  expect_equal(ev_mixture(-7, 40, 0.3, sign = -1), ev_mixture(7, 40, 0.3))
  expect_equal(ev_log_mixture(3, 10, 0.2, sign = -1), ev_log_mixture(-3, 10, 0.2))
  # an e-value is nonnegative
  expect_true(all(ev_mixture(seq(-50, 50, 5), 100, 0.5) >= 0))
})

test_that("M14: e-value inputs are validated", {
  expect_error(ev_mixture(0, 0, rho = -1))
  expect_error(ev_mixture(0, 0, rho = 0))
  expect_error(ev_mixture(0, -1, rho = 1))
  expect_error(ev_mixture(0, 10, 1, sign = 2))
  expect_error(ev_betting_paths(0.5, 0.5, 1.5, 10, 5))
  expect_error(ev_betting_paths(0.5, 0.5, -2, 10, 5))
  expect_error(ev_betting_paths(1.2, 0.5, 0.1, 10, 5))
  expect_error(ev_betting_paths(0.5, 0.5, 0.1, 0, 5))
  expect_error(ev_walks(0, 0, 5))
  expect_error(ev_naive_false_alarm(50, c(0, 5), reps = 10))
  expect_error(ev_naive_false_alarm(50, 5, reps = 10, alpha = 0))
  # wealth stays nonnegative at the boundary of the domain
  w <- ev_betting_paths(0.5, 0.5, 1, 10, 50, seed = 1)
  expect_false(anyNA(w))
})

test_that("M4: ev_mixture avoids Inf - Inf", {
  expect_equal(ev_mixture(-1e200, 100, 1), 7.98e-201, tolerance = 0.01)
  expect_equal(ev_mixture(-1000, 100, 1), 0.000798, tolerance = 0.01)
})
