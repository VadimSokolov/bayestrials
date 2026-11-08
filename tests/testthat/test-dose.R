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
  chk(dose_boin, "A", 0.551, 0.2945, 29.99)
  chk(dose_boin, "B", 0.4165, 0.491, 30)
})
