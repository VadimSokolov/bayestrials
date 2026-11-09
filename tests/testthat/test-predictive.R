test_that("normal-scale functions match Python", {
  expect_equal(cond_power(1.5, 0.5, 2), 0.5566203042214791, tolerance = 1e-6)
  expect_equal(pp_normal(1.5, 0.5), 0.5640936316683145, tolerance = 1e-6)
  expect_equal(pp_normal(1.5, 0.5, 0.5, 4), 0.45253698124206054, tolerance = 1e-6)
})

test_that("exact beta-binomial predictive probabilities match Python", {
  expect_equal(pp_one_arm(8, 20, 40, 0.2, 0.95), 0.8922090244562376, tolerance = 1e-6)
  expect_equal(pp_two_arm(8, 20, 5, 20, 10, 10, 0.9), 0.48164917532567814, tolerance = 1e-6)
})

test_that("Lee-Liu table and OC match Python", {
  tb <- lee_liu_table(20, 0.2, c(10, 15), 0.05, 0.95, 0.9)
  expect_equal(tb$s_fut, c(1, 3)); expect_equal(tb$s_eff, c(5, 7))
  oc <- lee_liu_oc(20, 0.2, tb, 0.2, 0.9)
  expect_equal(unname(oc), c(0.09410525053756386, 14.438103792025617), tolerance = 1e-6)
})

test_that("backward induction and forward OC match Python", {
  bi <- backward_induction(12, 0.01)
  expect_equal(sum(vapply(bi$stop, sum, 0)), 753)
  expect_equal(sum(vapply(bi$action, sum, 0)), 364)
  expect_equal(unname(forward_oc(bi, 0.5, 0.3)), c(0.6993814994591284, 1.9982328555975097),
               tolerance = 1e-6)
  bc <- backward_induction(12, 0.0005, "class", 0.9)
  expect_equal(sum(vapply(bc$stop, sum, 0)), 604)
  # 220, not the earlier grid value 221: at k = 7, (s1, s0) = (2, 0) the exact posterior
  # probability is exactly 0.9 = gamma, and the rule `q > gamma` is strict (the 600-node grid
  # put it at 0.9 + 1e-6)
  expect_equal(sum(vapply(bc$action, sum, 0)), 220)
  expect_equal(unname(forward_oc(bc, 0.5, 0.3)), c(0.3294893876600838, 8.804530390009772),
               tolerance = 1e-6)
})

test_that("fixed_oc matches Python", {
  expect_equal(fixed_oc(20, 0.95, 0.5, 0.3), 0.3243982367746252, tolerance = 1e-6)
})

test_that("M3: pp_two_arm and the state grid use exact superiority", {
  expect_equal(pp_two_arm(0, 0, 0, 0, 0, 0, eta = .1, pri = rep(.01, 4)), 1)
  expect_equal(pp_two_arm(0, 0, 0, 0, 0, 0, eta = .6, pri = rep(.01, 4)), 0)
  expect_equal(pp_two_arm(0, 0, 0, 0, 0, 0, eta = .4), 1)
  Q <- bayestrials:::.pr_superior_grid(.01 + 0:3, .01 + 3:0, .01 + 0:3, .01 + 3:0)
  expect_equal(Q + t(Q), matrix(1, 4, 4), tolerance = 1e-9)
  expect_equal(diag(Q), rep(.5, 4))
  Qn <- bayestrials:::.pr_superior_grid(c(3, 8), c(9, 4), c(2, 5, 10), c(6, 6, 2))
  expect_equal(Qn[2, 3], prob_superior(8, 4, 10, 2), tolerance = 1e-9)
  expect_equal(Qn[1, 2], prob_superior(3, 9, 5, 6), tolerance = 1e-9)
  # with no data the exchangeable posterior probability is exactly 1/2
  expect_equal(fixed_oc(0, 0.4, 0.3, 0.3), 1)
})
