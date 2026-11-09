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

test_that("capped allocation respects the floor in every history", {
  set.seed(3)
  R <- 400
  for (K in c(2, 3, 4)) {
    S <- matrix(rpois(R * K, 8), R, K) * matrix(rbinom(R * K, 1, 0.5) * 4 + 1, R, K)
    F <- matrix(rpois(R * K, 8), R, K)
    for (fl in c(0.1, 0.2)) {
      lower <- if (K > 2) fl / K else fl
      p <- rar_alloc_probs("capped", S, F, 200, 400, K, floor = fl)
      expect_true(all(p >= lower - 1e-9))
      expect_equal(rowSums(p), rep(1, R), tolerance = 1e-9)
      if (K == 2) expect_true(all(p <= 1 - fl + 1e-9))
      if (K > 2) {
        q <- rar_alloc_probs("capped", S, F, 200, 400, K, floor = fl, fixed_ctrl = TRUE)
        expect_equal(q[, 1], rep(1 / K, R))
        expect_true(all(q[, -1] >= lower - 1e-9))
        expect_equal(rowSums(q), rep(1, R), tolerance = 1e-9)
      }
    }
  }
  # extreme history: raw Thompson probability is (0, 1)
  p <- rar_alloc_probs("capped", matrix(c(0, 500), 1), matrix(c(500, 0), 1), 200, 400, 2, floor = 0.1)
  expect_equal(as.vector(p), c(0.1, 0.9))
  expect_equal(rar_bound_probs(matrix(c(0, 1), 1), 0.1), matrix(c(0.1, 0.9), 1))
  expect_equal(rar_bound_probs(matrix(c(0, 0.2, 0.8), 1), 0.1), matrix(c(0.1, 0.18, 0.72), 1), tolerance = 1e-6)
  expect_error(rar_alloc_probs("capped", matrix(0, 1, 2), matrix(0, 1, 2), 200, 400, 2, floor = 0.6))
})

test_that("fixed-control allocation uses probabilities of being best among experimental arms", {
  set.seed(4)
  # control is by far the best arm, experimental arms 2 and 3 are equal and arm 4 is worse
  S <- matrix(c(90, 20, 20, 5), 1); F <- matrix(c(10, 80, 80, 95), 1)
  q <- rar_alloc_probs("capped", S, F, 200, 400, 4, floor = 0.1, fixed_ctrl = TRUE, M = 20000)
  expect_equal(q[1], 0.25)
  expect_equal(q[2], q[3], tolerance = 0.02)
  expect_gt(q[2] / q[4], 10)
  expect_equal(sum(q[2:4]), 0.75)
  # with the control excluded the experimental arms share the full 0.75 (all-arm normalisation
  # would leave them with only the residual best-arm mass)
  S2 <- matrix(c(90, 50, 5, 5), 1); F2 <- matrix(c(10, 50, 95, 95), 1)
  q2 <- rar_alloc_probs("capped", S2, F2, 200, 400, 4, floor = 0.1, fixed_ctrl = TRUE, M = 20000)
  expect_gt(q2[2], 0.69)
})

test_that("rar_simulate enrols exactly N patients", {
  for (cfg in list(c(10, 20), c(101, 20), c(7, 3), c(100, 20), c(45, 50))) {
    for (rule in c("equal", "ts", "capped")) {
      r <- rar_simulate(c(1, 1), cfg[1], rule, reps = 30, batch = cfg[2], seed = 1)
      expect_equal(r$resp, cfg[1])
      r4 <- rar_simulate(rep(1, 4), cfg[1], rule, reps = 10, batch = cfg[2], seed = 1)
      expect_equal(r4$resp, cfg[1])
    }
  }
  r <- rar_simulate(c(1, 1), 101, "equal", reps = 50, batch = 20, seed = 2)
  expect_equal(r$n_inferior, 101 * (1 - r$share_best))
})

test_that("rar_simulate keeps allocation within bounds in every simulated batch", {
  seen <- list()
  orig <- rar_alloc_probs
  testthat::local_mocked_bindings(rar_alloc_probs = function(...) {
    p <- orig(...)
    seen[[length(seen) + 1L]] <<- p
    p
  })
  rar_simulate(c(0.1, 0.9), 203, "capped", reps = 200, batch = 20, seed = 5, floor = 0.1)
  expect_gt(length(seen), 5)
  expect_true(all(vapply(seen, function(p) all(p >= 0.1 - 1e-9 & p <= 0.9 + 1e-9), logical(1))))
  seen <- list()
  rar_simulate(c(0.1, 0.1, 0.9, 0.9), 203, "capped", reps = 100, batch = 20, seed = 5,
               fixed_ctrl = TRUE)
  expect_true(all(vapply(seen, function(p) all(p[, -1] >= 0.1 / 4 - 1e-9) && all(abs(p[, 1] - 0.25) < 1e-12), logical(1))))
})

test_that("best-arm rejection summaries are sane under the null and with a single replicate", {
  r <- rar_simulate(c(0.3, 0.3), 100, "equal", reps = 200, seed = 1)
  expect_false(anyNA(unlist(r)))
  expect_true(r$rej_best_naive >= 0 && r$rej_best_naive <= 0.1)
  r4 <- rar_simulate(rep(0.3, 4), 120, "ts", drift = 0.2, reps = 100, seed = 1)
  expect_false(anyNA(unlist(r4)))
  expect_true(r4$rej_best_strat <= r4$rej_any_strat + 1e-12)
  # control strictly best: summary still defined
  rc <- rar_simulate(c(0.6, 0.3, 0.3), 100, "equal", reps = 100, seed = 1)
  expect_false(anyNA(unlist(rc)))
  # alternative unchanged: best arm is the best experimental arm
  ra <- rar_simulate(c(0.3, 0.3, 0.4, 0.5), 200, "equal", reps = 300, seed = 1)
  expect_lte(ra$rej_best_naive, ra$rej_any_naive + 1e-12)
  d <- rar_drift_sweep(drifts = 0.1, reps = 20, N = 100, seed = 1)
  expect_false(anyNA(d$rej_best_naive))
  expect_false(anyNA(rar_four_arm(reps = 20, N = 100, seed = 1)$null$rej_best_naive))
  for (rule in c("equal", "ts", "capped")) {
    one <- rar_simulate(c(0.3, 0.45), 60, rule, reps = 1, seed = 1)
    expect_false(anyNA(unlist(one)))
    one4 <- rar_simulate(c(0.3, 0.3, 0.4, 0.5), 60, rule, reps = 1, seed = 1, fixed_ctrl = rule == "capped")
    expect_false(anyNA(unlist(one4)))
  }
  expect_equal(dim(rar_alloc_probs("ts", matrix(c(3, 6), 1), matrix(c(7, 4), 1), 20, 100, 2)), c(1, 2))
})

test_that("N2/M5: rar_simulate validates N and batch, floor feasibility up front", {
  expect_error(rar_simulate(c(.3, .4), 10, "equal", reps = 2, batch = 0))
  expect_error(rar_simulate(c(.3, .4), 10.5, "equal", reps = 2))
  expect_error(rar_simulate(c(.3, .4), 10, "equal", reps = 2, batch = 2.5))
  expect_error(rar_simulate(c(.3, .4), Inf, "equal", reps = 2))
  expect_error(rar_alloc_probs("capped", matrix(0, 1, 2), matrix(0, 1, 2), 0, 100, 2, floor = 0.6))
})
