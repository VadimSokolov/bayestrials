if (!exists("select_design")) load(testthat::test_path("..", "..", "data", "select_design.rda"), envir = globalenv())

test_that("look table matches Python and is benefit-positive", {
  lt <- select_look_table(952, 0.8)
  expect_equal(lt$events, c(238L, 476L, 714L, 952L))
  expect_equal(lt$z, c(1.721245522, 2.434208762, 2.981284697, 3.442491044), tolerance = 1e-6)
  expect_equal(lt$log_bf[4], 4.169747724, tolerance = 1e-6)
  expect_equal(lt$post_prob[1], 0.9367456861, tolerance = 1e-6)
  # harmful hazard ratio: signed statistic is negative and the probability of benefit is small
  lh <- select_look_table(952, 1 / 0.8)
  expect_equal(lh$z, -lt$z, tolerance = 1e-10)
  expect_equal(lh$post_prob, 1 - lt$post_prob, tolerance = 1e-10)
  expect_true(all(lh$post_prob < 0.5))
  expect_equal(select_post_prob(10, 238) + select_post_prob(-10, 238), 1)
})

test_that("month of events is monotone and hits endpoints", {
  expect_equal(select_month_of_events(1270, 1270), 60, tolerance = 1e-6)
  m <- select_month_of_events(c(100, 464, 900), 1270)
  expect_true(all(diff(m) > 0))
})

test_that("small SELECT simulation is in loose bounds", {
  r <- select_simulate(1270, 0.8, nsim = 1500, seed = 7, batch = 750)
  expect_gt(r$pstar, 0.985); expect_lt(r$pstar, 0.9995)
  expect_gt(r$naive_typeI, 0.07); expect_lt(r$naive_typeI, 0.15)
  expect_lt(abs(r$bayes_typeI - 0.025), 0.002)
  med <- r$rows$median_events_success
  expect_true(med[3] < med[2] && med[2] < med[1])
  expect_equal(med[3], 464, tolerance = 0.15)
  expect_equal(med[2], 608, tolerance = 0.15)
  expect_true(med[1] >= 600 && med[1] <= 960)
  expect_equal(r$gs_bounds, c(4.0487, 2.8628, 2.3375, 2.0243), tolerance = 2e-3)
})

test_that("group sequential comparators respect alpha and odd event totals", {
  r <- select_simulate(401, 0.8, nsim = 6000, alpha = 0.10, seed = 3, batch = 3000)
  expect_equal(r$dg[length(r$dg)], 401)
  expect_equal(r$gs_bounds, gs_boundary(4, 0.10, t = c(100, 200, 300, 401) / 401), tolerance = 0.05)
  expect_lt(abs(r$gs_typeI - 0.10), 0.02)
  expect_lt(abs(r$actual_typeI - 0.10), 0.02)
  expect_lte(r$ev_typeI, 0.10 + 0.01)  # Ville bound: conservative
  r0 <- select_simulate(401, 0.8, nsim = 6000, alpha = 0.025, seed = 3, batch = 3000)
  expect_lt(r0$gs_bounds[4], r$gs_bounds[4] + 1)  # sanity: bounds differ with alpha
  expect_true(all(r$gs_bounds < r0$gs_bounds))
  expect_lt(abs(r0$gs_typeI - 0.025), 0.01)
})

test_that("unconditional summaries include non-crossing trials", {
  r <- select_simulate(400, 0.9, nsim = 3000, seed = 5, batch = 1500)$rows
  expect_true(all(c("median_events_success", "mean_events_success", "median_events_all",
                    "mean_events_all", "power", "prob_cross") %in% names(r)))
  expect_true(all(r$power < 1 & r$power > 0))
  expect_true(all(r$median_events_all >= r$median_events_success))
  expect_true(all(r$mean_events_all >= r$mean_events_success))
  expect_equal(r$prob_cross, r$power)
  mx <- c(400, 400, 400, 400)
  # mean over all trials = p * conditional mean + (1 - p) * maximum
  expect_equal(r$mean_events_all, r$power * r$mean_events_success + (1 - r$power) * mx, tolerance = 1e-8)
  expect_true(all(r$mean_events_all <= mx))
})

test_that("actual SELECT design comparator", {
  expect_equal(select_design$design_events, 1225L)
  expect_equal(select_design$design_HR, 0.83)
  expect_equal(select_design$design_interim_events, 817L)
  expect_equal(select_design$observed_HR, 0.80)
  expect_equal(select_design$observed_events, 1270L)
  expect_equal(select_design$design_obf_bounds, gs_boundary(2, 0.025, t = c(817 / 1225, 1)), tolerance = 1e-4)
  r <- select_simulate(1270, 0.80, nsim = 3000, seed = 9, batch = 1500)
  expect_equal(r$actual_bounds, gs_boundary(2, 0.025, t = c(817 / 1225, 1)), tolerance = 1e-8)
  expect_equal(r$max_events[4], 1225)
  a <- r$rows[4, ]
  expect_match(a$rule, "SELECT actual design")
  expect_true(a$median_events_success %in% c(817, 1225))
  expect_lte(a$mean_events_all, 1225)
  # true HR 0.80 is better than the design HR 0.83: power above 90%
  expect_gt(a$power, 0.90)
  # only the two analyses can end the trial early
  expect_lt(abs(r$actual_typeI - 0.025), 0.01)
  # under the design HR the power is close to 90%
  r83 <- select_simulate(1270, 0.83, nsim = 4000, seed = 9, batch = 2000)$rows[4, ]
  expect_equal(r83$power, 0.90, tolerance = 0.05)
})
