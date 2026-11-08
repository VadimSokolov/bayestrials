test_that("rar_simulate agrees with Python within 4 MC se (two-arm)", {
  reps <- 3000
  se <- function(p) 4 * sqrt(2 * p * (1 - p) / reps)
  e <- rar_simulate(c(0.3, 0.45), 400, "equal", reps = reps, seed = 11)
  expect_lt(abs(e$rej_best_naive - 0.88), se(0.88))
  expect_lt(abs(e$resp - 149.77), 4 * sqrt(2) * 8.6 / sqrt(reps) * 1.5)
  t <- rar_simulate(c(0.3, 0.45), 400, "ts", reps = reps, seed = 11)
  expect_lt(abs(t$share_best - 0.8666), 0.03)
  expect_lt(abs(t$rej_best_naive - 0.6197), se(0.62))
  expect_lt(abs(t$rej_best_strat - 0.494), se(0.49))
  cp <- rar_simulate(c(0.3, 0.45), 400, "capped", reps = reps, seed = 11)
  expect_lt(abs(cp$share_best - 0.7844), 0.03)
  expect_lt(abs(cp$rej_best_naive - 0.7373), se(0.74))
  expect_lt(abs(cp$sd_n_best - 20.46), 3)
})

test_that("rar_simulate agrees with Python (four-arm null with drift)", {
  reps <- 1500
  se <- function(p) 4 * sqrt(2 * p * (1 - p) / reps)
  r <- rar_simulate(rep(0.3, 4), 600, "ts", drift = 0.2, reps = reps, seed = 14)
  expect_lt(abs(r$rej_any_naive - 0.1713), se(0.17))
  expect_lt(abs(r$rej_any_strat - 0.012), se(0.012) + 0.01)
  f <- rar_simulate(rep(0.3, 4), 600, "capped", drift = 0.2, reps = reps, seed = 14, fixed_ctrl = TRUE)
  expect_lt(abs(f$share_ctrl - 0.2504), 0.01)
  a <- rar_simulate(c(0.3, 0.3, 0.4, 0.5), 600, "capped", reps = reps, seed = 13, fixed_ctrl = TRUE)
  expect_lt(abs(a$share_best - 0.505), 0.03)
  expect_lt(abs(a$rej_best_naive - 0.932), se(0.93) + 0.01)
})

test_that("allocation probabilities are valid", {
  S <- matrix(c(5, 10, 2, 1), 2); F <- matrix(c(15, 10, 8, 9), 2)
  for (rule in c("equal", "ts", "capped")) {
    p <- rar_alloc_probs(rule, S, F, 200, 400, 2)
    expect_equal(rowSums(p), c(1, 1), tolerance = 1e-8)
  }
  p <- rar_alloc_probs("capped", cbind(S, S), cbind(F, F), 200, 400, 4, fixed_ctrl = TRUE)
  expect_equal(p[, 1], c(0.25, 0.25))
  expect_equal(rowSums(p), c(1, 1), tolerance = 1e-8)
})

test_that("wrappers return data.frames", {
  d <- rar_two_arm(reps = 20, N = 100, seed = 1)
  expect_equal(nrow(d), 3)
  expect_equal(nrow(rar_drift_sweep(drifts = c(0, 0.1), reps = 20, N = 100)), 6)
  expect_equal(nrow(rar_four_arm(reps = 20, N = 100)$alt), 4)
})

test_that("rar_gittins_index reproduces the chapter values", {
  expect_equal(round(rar_gittins_index(1, 1), 3), 0.703)
  expect_equal(round(rar_gittins_index(30, 70), 3), 0.308)
})

test_that("rar_regret_sim reproduces the chapter figures", {
  reg <- rar_regret_sim(T = 1000, reps = 300, seed = 3)
  expect_equal(dim(reg), c(1000, 3))
  expect_equal(colnames(reg), c("Equal", "UCB1", "Thompson"))
  expect_equal(unname(reg[100, "Equal"]), 50, tolerance = 0.1)
  expect_lt(abs(reg[100, "Thompson"] - 27), 4)
  expect_lt(reg[1000, "Thompson"], reg[1000, "UCB1"])
  expect_lt(reg[1000, "UCB1"], reg[1000, "Equal"])
})

test_that("rar_ecmo_urn reproduces the chapter numbers", {
  r <- rar_ecmo_urn(seed = 14)
  expect_lt(abs(r$n_cmt - 3.7), 0.1)
  expect_lt(abs(r$p_le1 - 0.11), 0.02)
  expect_lt(abs(r$power_rpw - 0.15), 4 * sqrt(0.15 * 0.85 / 2000))
  expect_lt(abs(r$power_fixed - 0.29), 4 * sqrt(0.29 * 0.71 / 2000))
})
