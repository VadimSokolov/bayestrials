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
