test_that("Lindley table matches zip type1_bayes_risk.tex:~110", {
  t <- decide_lindley_table()
  expect_equal(round(t$B, 3)[1:3], c(1.728, 0.666, 0.215), tolerance = 1e-3)
  expect_equal(round(t$B[4:5], 3), c(0.068, 0.022))
  expect_equal(round(t$post_null, 2), c(0.37, 0.60, 0.82, 0.94, 0.98))
})

test_that("exact Bayes bar carries (1+1/q); zip log(1+n) bar is approximate", {
  expect_equal(decide_bayes_bar(10), sqrt(1.1 * log(11)))
  expect_equal(decide_crossing_q(), 41.58, tolerance = 1e-3)
  t <- decide_bayes_type1_table()
  expect_equal(round(t$bar_approx, 2), c(1.55, 2.15, 2.63, 3.03))
  expect_equal(round(t$alpha_approx[2:4], 4), c(0.0317, 0.0085, 0.0024), tolerance = 0.02)
  expect_equal(round(t$alpha_exact, 4), c(0.1044, 0.0309, 0.0085, 0.0024))
  # exact bar at B = 1 is where decide_bf_q = 1
  expect_equal(decide_bf_q(decide_bayes_bar(100), 100), 1)
})

test_that("implied c, Bayes risk and loss ratios", {
  # the test |Z|>b is the Bayes rule at its implied c
  q <- 50; b <- 1.96; cc <- decide_implied_c(b, q)
  expect_equal(decide_bf_q(b, q), cc)
  r <- decide_bayes_risk(c(1.96, decide_bayes_bar(q)), q)
  expect_lt(r$risk[2], r$risk[1])
  expect_equal(decide_loss_ratio(c(0.975, 0.995, 0.99))$ratio, c(39, 199, 99), tolerance = 1e-9)
  expect_equal(decide_loss_ratio(ratio = 10)$pstar, 10 / 11)
})

test_that("Jeffreys table matches reg_evalues.tex:~63-80", {
  t <- decide_jeffreys_table()
  expect_equal(round(t$t_crit, 2), c(1.66, 2.25, 2.71, 3.11, 3.28, 3.46, 3.78))
  expect_equal(round(t$p_two_sided, 4)[1:3], c(0.0970, 0.0245, 0.0067), tolerance = 1e-2)
  expect_equal(round(t$t_strong, 2), c(2.71, 3.11, 3.46, 3.78, 3.92, 4.07, 4.34))
  expect_equal(round(t$K_exact_at_tcrit[c(1, 2, 3)], 3), c(1.302, 1.059, 1.008))
  expect_equal(round(decide_jeffreys_k_exact(sqrt(log(pi * 1e5 / 2)), 1e5), 4), 1.0001)
})

test_that("square-root rule matches reg_planner.tex:~60", {
  t <- size_sqrt_table(c(1e5, 1e6), 1.3, 2)
  expect_equal(t$n_star[2], 1088, tolerance = 1e-3)
  expect_equal(t$n_opt[2], 1084, tolerance = 2e-3)
  expect_equal(round(t$loss_min[2]), 2251, tolerance = 1e-3)
  expect_equal(round(t$loss_approx[2]), 2257, tolerance = 1e-3)
  expect_equal(t$rel_error[2], 0.003, tolerance = 0.3)
  expect_equal(size_sqrt_rule(1e6, 0.5 * 1, 1)[["n_star"]], 1414.2, tolerance = 1e-4)
  # cost remark
  expect_equal(size_sqrt_rule(1e6, 1, 1, cost = 0.1)[["n_star"]], sqrt(1e6 / (2 + 2 * sqrt(2 * pi) * 0.1)))
})

test_that("cost of K subsets matches reg_subsets.tex:29-34", {
  t <- decide_subset_cost()
  expect_equal(round(t$z_K, 2), c(1.96, 2.57, 2.80, 3.29, 3.48, 4.05))
  expect_equal(round(t$ratio, 2), c(1.00, 1.48, 1.69, 2.17, 2.38, 3.05))
  expect_equal(round(t$extra_per_1000), c(0, 485, 693, 1172, 1377, 2052), tolerance = 2e-3)
})

test_that("K = 20 simulation reproduces size and the negative power result", {
  s <- decide_max_of_k_sim(nsim = 1e5, seed = 7, power_rho = 0)
  expect_equal(s$size$naive[1], 0.398, tolerance = 0.03)
  expect_equal(s$size$sidak[1], 0.025, tolerance = 0.1)
  expect_equal(s$size$naive[2:3], c(0.210, 0.112), tolerance = 0.06)
  expect_lt(s$size$mixture[1], 0.002)
  expect_equal(s$power$sidak[1:2], c(0.505, 0.692), tolerance = 0.05)
  expect_true(all(s$power$mixture < s$power$sidak))
})

test_that("direction = 'positive' reproduces the chapter's mk_dir (11-evalues.qmd)", {
  mk_dir <- function(K = 20, rho = c(0, 0.5, 0.8), alpha = 0.025, lam2 = 9, nsim = 2e5,
                     scenarios = list(`one effect d=3.0` = 3, `three effects d=2.5` = rep(2.5, 3)),
                     power_rho = c(0, 0.5), seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    c0 <- (1 + lam2)^-0.5; bb <- lam2 / (2 * (1 + lam2))
    draw <- function(theta, r) sweep(sqrt(r) * rnorm(nsim) + sqrt(1 - r) * matrix(rnorm(nsim * K), nsim, K), 2, theta, "+")
    stat <- function(z) list(mx = do.call(pmax, as.data.frame(z)), E = rowMeans(2 * c0 * exp(bb * z^2) * (z > 0)))
    zero <- numeric(K); size <- list(); cal <- list()
    for (r in rho) {
      s <- stat(draw(zero, r))
      size[[length(size) + 1]] <- data.frame(rho = r, naive = mean(s$mx > qnorm(1 - alpha)), sidak = mean(s$mx > qnorm((1 - alpha)^(1 / K))), mixture = mean(s$E >= 1 / alpha))
      cal[[as.character(r)]] <- c(sidak = unname(quantile(s$mx, 1 - alpha)), mixture = unname(quantile(s$E, 1 - alpha)))
    }
    pw <- list()
    for (r in power_rho) {
      th <- cal[[as.character(r)]]
      for (nm in names(scenarios)) {
        t <- zero; t[seq_along(scenarios[[nm]])] <- scenarios[[nm]]
        s <- stat(draw(t, r))
        pw[[length(pw) + 1]] <- data.frame(scenario = nm, rho = r, sidak = mean(s$mx > th[["sidak"]]), mixture = mean(s$E >= th[["mixture"]]))
      }
    }
    list(size = do.call(rbind, size), power = do.call(rbind, pw))
  }
  sc <- list(`one effect d=3.0` = 3, `three effects d=2.5` = rep(2.5, 3))
  a <- mk_dir(K = 20, nsim = 2e4, scenarios = sc, seed = 7)
  b <- decide_max_of_k_sim(K = 20, nsim = 2e4, scenarios = sc, seed = 7, direction = "positive")
  expect_equal(b$size, a$size)
  expect_equal(b$power, a$power)
  d <- decide_max_of_k_sim(K = 20, nsim = 2e4, scenarios = sc, seed = 7)
  expect_false(isTRUE(all.equal(d$size$mixture, b$size$mixture)))
})

test_that("sample size criteria match decision_foundations.tex:~98-108", {
  expect_equal(decide_ss_alc(0.3), 342)
  p <- decide_ss_posterior_prob(c(200, 300, 400))$prob
  expect_equal(round(p, 2), c(0.77, 0.86, 0.91))
  p2 <- decide_ss_posterior_prob(c(50, 400), analysis_mean = 0.3, analysis_sd = 0.1)$prob
  expect_true(all(p2 > 0.985))
  expect_equal(round(decide_ss_posterior_prob(234, design_sd = 0)$prob, 2), 0.90, tolerance = 0.01)
  expect_equal(round(decide_ss_posterior_prob(234)$prob, 2), 0.81)
  # zip says 230; the exact crossing is 226 (zip is rounded)
  expect_equal(decide_ss_min_n(0.8), 226)
  a <- decide_ss_acc_woc(342, 0.3, design_sd = 0.1)
  expect_gt(a$ACC, 0)
  m <- decide_ss_acc_woc(100, 0.4, analysis_mean = 0.3, analysis_sd = 0.1)
  expect_equal(m$ACC, m$WOC, tolerance = 1e-6)
})

test_that("stopping table matches decision_foundations.tex:~125", {
  r <- decide_stop_single_arm(p = c(0.2, 0.3, 0.45))
  expect_equal(r$table$cont_lo, c(0, 1, 4, 7, 10, 14, 15))
  expect_equal(r$table$cont_hi, c(4, 5, 8, 11, 14, 16, 17))
  expect_equal(r$table$adopt_min, c(2, 4, 7, 10, 13, 16, 17))
  expect_equal(r$oc$P_adopt[1], 0.05, tolerance = 0.5)
  expect_equal(r$oc$P_adopt[3], 0.98, tolerance = 0.02)
  expect_equal(r$oc$E_n[3], 19, tolerance = 0.05)
  expect_equal(r$oc$E_n[2], 32, tolerance = 0.06)
  expect_equal(r$oc$E_n[1], 25, tolerance = 0.06)
})

test_that("pivotal sizing matches programme.tex:~58-66", {
  f <- decide_pivotal_table(m = 0.30, v = 0.141^2)
  expect_equal(round(f$assurance, 2), c(0.73, 0.80, 0.85, 0.87, 0.89, 0.91, 0.92))
  expect_equal(round(f$enpv), c(593, 644, 662, 665, 660, 640, 612), tolerance = 2e-3)
  expect_equal(round(f$power, 2)[1:2], c(0.85, 0.96))
  s <- decide_pivotal_table(m = 0.20, v = 0.115^2)
  expect_equal(round(s$assurance, 2), c(0.51, 0.61, 0.68, 0.72, 0.75, 0.79, 0.82))
  expect_equal(round(s$enpv), c(399, 470, 507, 526, 535, 532, 516), tolerance = 3e-3)
  expect_equal(round(s$power, 2)[1:3], c(0.52, 0.69, 0.81))
})

test_that("interim posterior and payer example match programme.tex", {
  r <- decide_interim_posterior(0.20, 0.115, 0.15, 0.082)
  expect_equal(c(r$m, r$s), c(0.167, 0.067), tolerance = 5e-3)
  expect_equal(r$pp, 0.79, tolerance = 0.01)
  i <- decide_inb(0.22, 0.058, 0.25, 1e5, 4000)
  expect_equal(c(i$mean, i$sd), c(1500, 1450))
  expect_equal(i$prob_positive, 0.85, tolerance = 0.01)
  expect_equal(decide_inb(0.22, 0.058, 0.25, 1e5, 3600)$prob_positive, 0.90, tolerance = 0.01)
})

test_that("NI table matches borrowing_examples.tex:~20-24", {
  t <- ni_posterior_table(0.01, 0.015, 0.04)
  expect_equal(round(t$post_mean, 3), c(0.010, 0.021, 0.003))
  expect_equal(round(t$prob_ni[1:2], 3), c(0.977, 0.945))
  expect_gt(t$prob_ni[3], 0.999)
  expect_equal(attr(t, "upper_limit"), 0.0394, tolerance = 1e-3)
  s <- ni_synthesis(0.01, 0.015, 0.10, 0.02, 0.5)
  expect_true(s$prob_retain > 0 && s$prob_retain < 1)
  expect_equal(s$fixed_margin, 0.5 * (0.10 - 1.96 * 0.02), tolerance = 1e-3)
})

test_that("rare disease and pediatric examples match borrowing_examples.tex", {
  t <- borrow_rare_table()
  expect_equal(round(t$p_superior, 2), c(0.83, 0.92, 0.95, 0.97))
  expect_equal(round(t$p_delta, 2), c(0.62, 0.74, 0.80, 0.86))
  expect_equal(t$ess_control, c(15, 25, 35, 55))
  expect_equal(attr(borrow_marginal_a0(yc = 6), "a0_hat"), 1)
  expect_equal(attr(borrow_marginal_a0(yc = 11), "a0_hat"), 0)
  p <- borrow_pediatric_mixture()
  expect_equal(c(p$weight, p$p_benefit, p$p_benefit_alone), c(0.87, 0.98, 0.88), tolerance = 0.01)
})
