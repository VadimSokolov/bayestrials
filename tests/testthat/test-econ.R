test_that("exposure after a decision follows the share of events", {
  x <- exposure_after_decision(635, 1270, 17604, 39.8 / 12, 701, 569)
  expect_equal(x$share, 0.5)
  expect_equal(x$placebo_py, 0.5 * 17604 * 39.8 / 12 / 2)
  expect_equal(x$placebo_events, 350.5)
  expect_equal(x$excess, 0.5 * (701 - 569))
  expect_equal(x$excess, x$placebo_events - 0.5 * 569)
  expect_equal(exposure_after_decision(1270, 1270, 100, 1, 10, 5)$excess, 0)
})

test_that("cost saving and label revenue are products of their factors", {
  expect_equal(trial_cost_saving(0.7, 157, 0.6), 0.3 * 0.6 * 157)
  expect_equal(trial_cost_saving(1, 157, 0.6), 0)
  expect_equal(label_revenue_forward(8, c(0.1, 0.25), 2), c(1.6, 4))
  expect_equal(label_revenue_forward(8, 0.1, 0), 0)
})
