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
