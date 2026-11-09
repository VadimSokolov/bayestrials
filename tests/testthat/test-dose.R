test_that("dose_boin_bounds matches Python and the BOIN package", {
  b <- dose_boin_bounds()
  expect_equal(unname(b), c(0.19680087055548712, 0.2983921523754597), tolerance = 1e-6)
  skip_if_not_installed("BOIN")
  for (tg in c(0.2, 0.25, 0.3)) {
    g <- BOIN::get.boundary(tg, 10, 3, p.saf = 0.6 * tg, p.tox = 1.4 * tg)
    b <- dose_boin_bounds(tg)
    expect_equal(unname(b[["escalate"]]), g$lambda_e, tolerance = 1e-6)
    expect_equal(unname(b[["deescalate"]]), g$lambda_d, tolerance = 1e-6)
  }
})

test_that("pava and true mtd match Python", {
  expect_equal(dose_pava(c(.1, .3, .2, .5), c(1, 1, 2, 1)), c(0.1, 0.7 / 3, 0.7 / 3, 0.5), tolerance = 1e-8)
  expect_equal(dose_true_mtd(dose_scenarios()$A), 3)
  expect_equal(dose_true_mtd(c(0.3, 0.5)), 0)
})

test_that("design simulations agree with Python within 4 MC se", {
  sc <- dose_scenarios()
  chk <- function(design, s, pc, pu, n_py) {
    reps <- 2000
    r <- dose_simulate(design, sc[[s]], reps, seed = 5)
    expect_lt(abs(r$p_correct - pc), 4 * sqrt(pc * (1 - pc) / reps) * sqrt(2))
    expect_lt(abs(r$p_under - pu), 4 * sqrt(pu * (1 - pu) / reps) * sqrt(2))
    expect_lt(abs(r$n - n_py), 1)
  }
  chk(dose_three_plus_three, "A", 0.372, 0.464, 13.887)
  chk(dose_three_plus_three, "B", 0.3565, 0.5225, 18.975)
  chk(dose_crm, "A", 0.663, 0.1485, 29.99)
  chk(dose_crm, "B", 0.528, 0.385, 30)
  # BOIN values re-derived after switching final selection to the BOIN::select.mtd rule
  # (they agree with BOIN::get.oc, see the next test); the previous Python values used count weights
  chk(dose_boin, "A", 0.5919, 0.2269, 29.99)
  chk(dose_boin, "B", 0.5136, 0.366, 30)
})

test_that("BOIN final selection agrees with BOIN::select.mtd", {
  skip_if_not_installed("BOIN")
  ref <- function(n, y, target = 0.25) {
    m <- BOIN::select.mtd(target, n, y)$MTD
    if (m == 99) 0L else as.integer(m)
  }
  hist <- list(
    list(n = c(6, 6, 3), y = c(0, 1, 0)),
    list(n = c(3, 3, 3), y = c(0, 0, 0)),
    list(n = c(3, 6, 6, 3), y = c(0, 1, 2, 2)),
    list(n = c(3, 3, 0, 0), y = c(0, 1, 0, 0)),
    list(n = c(6, 9, 6, 0, 0), y = c(1, 2, 3, 0, 0)),
    list(n = c(3, 3, 3), y = c(3, 0, 0)),
    list(n = c(9, 6, 3), y = c(3, 2, 3)),
    list(n = c(3, 3, 6, 12), y = c(0, 1, 2, 4)),
    list(n = c(6, 6, 6, 6), y = c(2, 1, 3, 2)))
  for (h in hist) expect_equal(dose_boin_select(h$n, h$y), ref(h$n, h$y))
  expect_equal(dose_boin_select(c(6, 6, 3), c(0, 1, 0)), 3L)
  expect_equal(dose_boin_select(c(3, 3, 3), c(0, 0, 0)), 3L)
  set.seed(7)
  for (i in 1:300) {
    k <- sample(3:6, 1)
    n <- 3 * rpois(k, 1.5); y <- rbinom(k, n, runif(1, 0.05, 0.6))
    if (sum(n) == 0) next
    tg <- sample(c(0.2, 0.25, 0.3), 1)
    expect_equal(dose_boin_select(n, y, tg), ref(n, y, tg))
  }
  # a trial with no DLTs selects the highest dose (tie rule below target)
  set.seed(1)
  expect_equal(dose_boin(rep(0, 3), N = 9)$selected, 3L)
})

test_that("3+3 and CRM invariants and skeleton validation", {
  set.seed(2)
  for (i in 1:50) {
    r <- dose_three_plus_three(c(0.05, 0.1, 0.25, 0.4))
    expect_true(all(r$n %% 3 == 0))
    expect_true(r$selected >= 0 && r$selected <= 4)
  }
  expect_equal(dose_three_plus_three(rep(0, 4))$selected, 4L)
  expect_equal(dose_three_plus_three(rep(1, 4))$selected, 0L)
  expect_error(dose_crm(c(0.1, 0.2, 0.3), N = 9), "skeleton")
  r <- dose_crm(c(0.05, 0.1, 0.25, 0.4, 0.55, 0.7), N = 10)
  expect_equal(sum(r$n), 9)
})

test_that("dose_boin operating characteristics agree with BOIN::get.oc", {
  skip_if_not_installed("BOIN")
  sc <- dose_scenarios()
  for (s in c("A", "B")) {
    r <- dose_simulate(dose_boin, sc[[s]], 4000, seed = 5)
    o <- BOIN::get.oc(0.25, sc[[s]], 10, 3, n.earlystop = 100, ntrial = 4000, seed = 1,
                      cutoff.eli = 0.95)
    expect_lt(max(abs(r$sel[-1] - o$selpercent / 100)), 0.04)
  }
})
