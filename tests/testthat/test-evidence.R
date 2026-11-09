test_that("evidence functions match Python to 1e-6", {
  expect_equal(false_positive_risk(.5, .8, .025), 0.030303030303030304, tolerance = 1e-6)
  expect_equal(power_one_sided(2.8), 0.7995559032981122, tolerance = 1e-6)
  expect_equal(bf10_normal(2.5, 10, log = TRUE), 1.6419614545099055, tolerance = 1e-6)
  expect_equal(posterior_null(1.3, .5), 0.2141650169574414, tolerance = 1e-6)
  expect_equal(min_bf(0.01), 0.1251815043353279, tolerance = 1e-6)
  expect_equal(min_bf(2, "gaussian"), 0.1353352832366127, tolerance = 1e-6)
  expect_equal(bayes_bar(10, 4), 5.687532394541968, tolerance = 1e-6)
  expect_equal(crossing_n(), 41.57329443857865, tolerance = 1e-6)
  expect_equal(crossing_n(0.05), 10.925408185346662, tolerance = 1e-6)
})

test_that("posterior-criterion calibration matches Python", {
  expect_equal(type1_posterior_flat(.1, .95, .05), 0.0001337720145987431, tolerance = 1e-6)
  expect_equal(type1_posterior_normal(.1, .95, .05, .3), 9.839658829624565e-05, tolerance = 1e-6)
  expect_equal(calibrate_c_flat(.1, .05), 0.4840322065576682, tolerance = 1e-6)
  expect_equal(calibrate_c_normal(.1, .05, .3), 0.46243902142783316, tolerance = 1e-6)
  expect_equal(power_posterior_normal(.1, .95, .05, .3, .15), 0.23480987260062558, tolerance = 1e-6)
  expect_equal(threshold_flat(.1, .95, .05), 0.1822426813475736, tolerance = 1e-6)
  expect_equal(type1_posterior_flat(0, calibrate_c_flat(0, .05, .025), .05), 0.025, tolerance = 1e-9)
})

test_that("M13: crossing_n finds both roots of a non-monotone Bayes bar", {
  r <- crossing_n(0.025, 2, all = TRUE)
  expect_equal(r, c(1.4894167, 5.3781852), tolerance = 1e-6)
  z2 <- qnorm(0.975)^2
  expect_equal(bayes_bar(r, 2), rep(z2, 2), tolerance = 1e-10)
  expect_equal(crossing_n(0.025, 2), r[2], tolerance = 1e-10)
  # monotone cases: one root, equation residual zero
  for (cc in c(0.5, 1)) {
    k <- crossing_n(0.025, cc)
    expect_equal(bayes_bar(k, cc), z2, tolerance = 1e-10)
    expect_length(crossing_n(0.025, cc, all = TRUE), 1)
  }
  # the dip stays above the fixed bar: no crossing
  expect_true(is.na(crossing_n(0.025, 50)))
  expect_length(crossing_n(0.025, 50, all = TRUE), 0)
  expect_error(crossing_n(alpha = 2))
  expect_error(crossing_n(c = -1))
})

test_that("M4: posterior_null is stable on extreme log Bayes factors", {
  expect_equal(posterior_null(-1000), 1)
  expect_equal(posterior_null(1000), 0)
  expect_equal(posterior_null(c(-1000, 0, 1000), 0.5), c(1, 0.5, 0))
  expect_equal(posterior_null(1.3, 0.2), 1 / (1 + 4 * exp(1.3)), tolerance = 1e-12)
  expect_equal(posterior_null(log(3), 0), 0)
  expect_equal(posterior_null(0, 1), 1)
  expect_error(posterior_null(0, 1.2))
})

test_that("MINOR 6: min_bf domain", {
  expect_equal(min_bf(0), 0)
  expect_equal(min_bf(1 / exp(1)), 1, tolerance = 1e-12)
  expect_equal(min_bf(c(0.5, 1)), c(1, 1))
  expect_equal(min_bf(c(0.001, 0.2)), -exp(1) * c(0.001, 0.2) * log(c(0.001, 0.2)))
  expect_true(all(diff(min_bf(seq(0, 0.9, length.out = 50))) >= 0))
  expect_true(is.na(min_bf(NA_real_)))
  expect_error(min_bf(-0.1))
  expect_error(min_bf(1.1))
  expect_equal(min_bf(0, "gaussian"), 1)
})

test_that("MINOR 7 and M14: false_positive_risk uses power and validates", {
  fpr <- function(pi, pw, a) a * (1 - pi) / (a * (1 - pi) + pw * pi)
  expect_equal(false_positive_risk(0.1, 0.9, 0.05), fpr(0.1, 0.9, 0.05))
  expect_equal(false_positive_risk(0.5, 1, 0.025), 0.025 / 1.025)
  expect_error(false_positive_risk(1.5, 0.8))
  expect_error(false_positive_risk(0.5, -0.1))
})
