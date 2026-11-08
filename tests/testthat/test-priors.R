test_that("prob_success and mixture match Python", {
  expect_equal(prob_success(0.3, 0.1, 0.2), 0.3403880475951664, tolerance = 1e-6)
  mp <- mixture_posterior(0.2, 0.1, data.frame(weight = c(.5, .5), mean = c(0, .3), sd = c(.5, .2)))
  expect_equal(mp$weights, c(0.30975768, 0.69024232), tolerance = 1e-6)
  expect_equal(mp$means, c(0.19230769, 0.22), tolerance = 1e-6)
  expect_equal(mp$sds, c(0.09805807, 0.08944272), tolerance = 1e-6)
  expect_equal(mp$p_pos, 0.9874783520863557, tolerance = 1e-6)
})
test_that("ess_events and cmp_samples", {
  expect_equal(ess_events(0.5), 16)
  d <- cmp_samples(1, 4, 4, 1, 10, 50, 2000, seed = 1)
  expect_equal(nrow(d), 2000)
  expect_true(mean(d$slope) > 0)
})
