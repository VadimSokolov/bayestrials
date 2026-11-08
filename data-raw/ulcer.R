# Run from the book root: Rscript bayestrials/data-raw/ulcer.R
ulcer <- data.frame(
  operation = c("Drainage and vagotomy", "25 per cent resection and vagotomy",
                "50 per cent resection and vagotomy", "75 per cent resection"),
  n = c(96L, 104L, 110L, 107L),
  recurrences = c(20L, 16L, 12L, 9L),
  stringsAsFactors = FALSE)
save(ulcer, file = "bayestrials/data/ulcer.rda", compress = "xz")

select_design <- list(
  alpha = 0.025, tau = 0.25, d_tot = 1270L, HR = 0.80,
  fractions = c(0.25, 0.5, 0.75, 1),
  obf_bounds = c(4.0487, 2.8628, 2.3375, 2.0243),
  info_per_event = 0.25, T_en = 36, T_end = 60,
  rewind_d_tot = 1257L, rewind_HR = 0.88)
save(select_design, file = "bayestrials/data/select_design.rda", compress = "xz")
