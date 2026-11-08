test_that("sim_confounded_trial matches Python within MC error", {
  set.seed(1)
  r <- t(replicate(4000, sim_confounded_trial(400)))
  expect_lt(abs(mean(r[, 1]) - 0.30715), 4 * 0.00075 * sqrt(2))
  expect_lt(abs(mean(r[, 2]) - 1.06137), 4 * 0.00134 * sqrt(2))
  r <- t(replicate(4000, sim_confounded_trial(400, assign = "coin")))
  expect_lt(abs(mean(r[, 1]) - 0.09745), 4 * 0.00078 * sqrt(2))
  expect_lt(abs(mean(r[, 2]) + 0.00173), 4 * 0.0016 * sqrt(2))
})

test_that("sim_adjusted_vs_naive matches Python", {
  set.seed(2)
  r <- t(replicate(4000, sim_adjusted_vs_naive(400)))
  expect_lt(abs(mean(r[, 1]) - 0.11853), 4 * 0.00095 * sqrt(2))
  expect_lt(abs(mean(r[, 2]) - 0.30656), 4 * 0.00074 * sqrt(2))
})

test_that("assignment functions behave", {
  a <- assign_block(100, sizes = 4, seed = 1)
  expect_length(a, 100)
  expect_true(all(tapply(a, rep(1:25, each = 4), sum) == 2))
  s <- rep(1:3, length.out = 60)
  as <- assign_stratified(s, seed = 1)
  expect_true(all(tapply(as, s, sum) == 10))
  expect_true(all(assign_simple(50, seed = 1) %in% 0:1))
  set.seed(3)
  X <- cbind(sample(0:1, 200, TRUE), sample(0:2, 200, TRUE))
  am <- assign_minimization(X, seed = 3)
  expect_lt(abs(sum(am) - 100), 10)
  expect_equal(assign_block(20, seed = 5), assign_block(20, seed = 5))
})

test_that("covariate_imbalance, arnold, compress match Python", {
  expect_equal(covariate_imbalance(c(0, 1, 0, 1), c(1, 1, 0, 0)), 0)
  a <- arnold_stochasticity(c(1, 3, 6, 7.5), 10)
  expect_equal(unname(a["S"]), 1.1, tolerance = 1e-6)
  expect_equal(unname(a["expected"]), 1.6, tolerance = 1e-6)
  expect_lt(abs(compress_ratio(rep(0:1, 50)) - 0.13), 0.03)
})

test_that("perm_test_arm is a valid p-value", {
  y <- c(3, 5, 8, 9, 12, 1, 2, 4, 6, 7); arm <- rep(1:0, each = 5)
  r <- perm_test_arm(y, arm, nsim = 20000, seed = 1)
  p_ex <- mean(apply(combn(10, 5), 2, function(i) mean(y[i]) - mean(y[-i])) >= r["obs"] - 1e-12)
  expect_lt(abs(r["p"] - p_ex), 4 * sqrt(p_ex * (1 - p_ex) / 20000))
  rb <- perm_test_arm(y, arm, nsim = 500, restricted_blocks = list(1:5, 6:10), seed = 1)
  expect_true(rb["p"] > 0 && rb["p"] <= 1)
})

test_that("rand_bits_lcg matches the direct recursion", {
  x <- 1; for (i in 2:50) x[i] <- (75 * x[i - 1]) %% 65537
  expect_equal(rand_bits_lcg(50), as.integer(x > 32768))
  expect_length(rand_bits_lcg(2000), 2000)
  expect_lt(compress_ratio(rand_bits_lcg(2000)), 0.25)
  expect_gt(compress_ratio(rand_bits_lcg(2000)), 0.125)
})
