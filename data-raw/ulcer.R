# Run from the book root: Rscript bayestrials/data-raw/ulcer.R
ulcer <- data.frame(
  operation = c("Drainage and vagotomy", "25 per cent resection and vagotomy",
                "50 per cent resection and vagotomy", "75 per cent resection"),
  n = c(96L, 104L, 110L, 107L),
  recurrences = c(20L, 16L, 12L, 9L),
  stringsAsFactors = FALSE)
save(ulcer, file = "bayestrials/data/ulcer.rda", compress = "xz")

# SELECT (Lincoff et al. 2023, NEJM 389, 2221-2232). Benefit-positive convention theta = -log(HR).
# Actual design: HR 0.83, 90% power, one-sided 0.025, 1225 planned first-event totals, one superiority
# interim at two thirds of the planned events, then the final analysis. Observed: HR 0.80, 1270 events.
# The first-event totals, HR and four-look constants below describe the book's hypothetical
# four-look reconstruction run at the observed 1270 events and HR 0.80.
design_interim <- as.integer(round(1225 * 2 / 3))
design_bounds <- c(2.4525, 2.0028)  # gs_boundary(2, 0.025, t = c(817 / 1225, 1)), O'Brien-Fleming type
select_design <- list(
  alpha = 0.025, tau = 0.25, d_tot = 1270L, HR = 0.80,
  fractions = c(0.25, 0.5, 0.75, 1),
  obf_bounds = c(4.0487, 2.8628, 2.3375, 2.0243),
  info_per_event = 0.25, T_en = 36, T_end = 60,
  rewind_d_tot = 1257L, rewind_HR = 0.88,
  design_HR = 0.83, design_power = 0.90, design_alpha = 0.025, design_events = 1225L,
  design_interim_fraction = 2 / 3, design_interim_events = design_interim,
  design_obf_bounds = design_bounds,
  observed_HR = 0.80, observed_events = 1270L)
save(select_design, file = "bayestrials/data/select_design.rda", compress = "xz")
