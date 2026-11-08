if (!exists("select_design")) load(testthat::test_path("..", "..", "data", "select_design.rda"), envir = globalenv())

test_that("look table matches Python", {
  lt <- select_look_table(952, 0.8)
  expect_equal(lt$events, c(238L, 476L, 714L, 952L))
  expect_equal(lt$z, c(1.721245522, 2.434208762, 2.981284697, 3.442491044), tolerance = 1e-6)
  expect_equal(lt$log_bf[4], 4.169747724, tolerance = 1e-6)
  expect_equal(lt$post_prob[1], 0.9367456861, tolerance = 1e-6)
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
  med <- r$rows$median_events
  expect_true(med[3] < med[2] && med[2] < med[1])
  expect_equal(med[3], 464, tolerance = 0.15)
  expect_equal(med[2], 608, tolerance = 0.15)
  expect_true(med[1] >= 600 && med[1] <= 960)
})
