test_that("core boundaries match published values", {
  expect_equal(gs_boundary(4), c(4.049, 2.863, 2.337, 2.024), tolerance = 1e-3)
  expect_equal(gs_boundary(4, type = "pocock"), rep(2.361, 4), tolerance = 1e-3)
  expect_equal(gs_spending_boundary((1:4) / 4), c(4.333, 2.963, 2.359, 2.014), tolerance = 1e-3)
})

test_that("two-sided recursion matches Python", {
  t <- (1:4) / 4
  b <- gs_boundary(4)
  expect_equal(gs_power(t, b), 0.025, tolerance = 1e-6)
  expect_equal(gs_power_drift(t, b), 3.2772396, tolerance = 1e-5)
  lo <- c(-9, 0, 0.5, 1.0)
  expect_equal(gs_power_drift(t, b, lo), 3.2898853, tolerance = 1e-5)
  expect_equal(gs_expected_information(t, b, lo, 0), 0.6920757, tolerance = 1e-5)
  expect_equal(gs_expected_information(t, b, lo, 3), 0.7782860, tolerance = 1e-5)
  expect_equal(gs_expected_information(t, b, NULL, 3), 0.7893430, tolerance = 1e-5)
  r <- gs_crossing_two_sided(t, b, lo, 3)
  expect_equal(r$up, c(0.00540795, 0.22401125, 0.37834727, 0.23227054), tolerance = 1e-6)
  expect_equal(r$down[2:4], c(1.69474267e-02, 1.03675129e-02, 1.09818068e-02), tolerance = 1e-6)
  r0 <- gs_crossing_two_sided(t, rep(1.96, 4))
  expect_equal(r0$up, c(0.0249979, 0.01655768, 0.01206932, 0.00945994), tolerance = 1e-6)
  expect_equal(gs_crossing_prob(rep(1.96, 4), t)$cross, r0$up, tolerance = 1e-5)
})

test_that("naive inflation and tolls", {
  expect_lt(abs(gs_naive_inflation(4) - 0.0638), 2e-3)
  expect_equal(gs_naive_inflation(1), 0.025, tolerance = 1e-3)
  expect_lt(abs(gs_naive_inflation(2) - 0.042), 2e-3)
  tl <- gs_tolls(c(25, 4000), tau = 0.25)
  expect_equal(tl$surcharge[1], 7.378, tolerance = 1e-4)
  expect_equal(tl$ockham, log(1 + c(25, 4000) * 0.0625))
  expect_equal(tl$anytime_bar^2 - tl$bayes_bar^2, rep(tl$surcharge[1], 2))
})
