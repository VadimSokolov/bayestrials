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
  expect_equal(sum(vapply(bc$action, sum, 0)), 221)
  expect_equal(unname(forward_oc(bc, 0.5, 0.3)), c(0.3294893876600838, 8.804530390009772),
               tolerance = 1e-6)
})

test_that("fixed_oc matches Python", {
  expect_equal(fixed_oc(20, 0.95, 0.5, 0.3), 0.3243982367746252, tolerance = 1e-6)
})
