#' Novick-Grizzle duodenal ulcer trial: four operations
#'
#' Recurrence counts for four operations for duodenal ulcer, used in the multiplicity chapter. The recurrence counts are illustrative: they are set at the scale of the 1965 trial and not transcribed from Novick and Grizzle; they are the counts in Polson and Sokolov's working paper, which itself states that they are illustrative.
#'
#' @format A data.frame with 4 rows and columns `operation`, `n` (patients), `recurrences`.
#' @source Illustrative counts after Novick and Grizzle (1965), Journal of the American Statistical Association 60, 81-96; Polson and Sokolov (2026), working paper.
"ulcer"

#' SELECT design constants
#'
#' Constants of the SELECT case study in the cases chapter, benefit-positive convention `theta = -log(HR)`. Two groups of fields. (1) The book's hypothetical reconstruction, run at the observed data: 1270 first-event totals, hazard ratio 0.80, prior sd 0.25 for `-log HR`, four equally spaced looks with one-sided 0.025 O'Brien-Fleming bounds on the Z scale, information 1/4 per event, linear enrolment over 36 months in a 60 month trial. (2) SELECT's actual design: hazard ratio 0.83 with 90% power at one-sided 0.025 and 1,225 planned events, one superiority interim at two thirds of the planned events (817) with O'Brien-Fleming-type bounds, then the final analysis; observed hazard ratio 0.80 and 1,270 events.
#'
#' @format A list with `alpha`, `tau`, `d_tot`, `HR`, `fractions`, `obf_bounds`, `info_per_event`, `T_en`, `T_end`, `rewind_d_tot`, `rewind_HR`, and the actual-design fields `design_HR`, `design_power`, `design_alpha`, `design_events`, `design_interim_fraction`, `design_interim_events`, `design_obf_bounds` (two-look bounds at information fractions 817/1225 and 1), `observed_HR`, `observed_events`.
#' @source SELECT: Lincoff et al. (2023), NEJM 389, 2221-2232 (semaglutide, observed HR 0.80); REWIND: Gerstein et al. (2019), Lancet 394, 121-130 (HR 0.88, 1257 events); bounds from Polson and Sokolov's calibration code.
"select_design"

#' ClinicalTrials.gov counts by year and term
#'
#' Counts of interventional studies starting each year (2000 to 2025) that mention Bayesian or adaptive terms, with the denominator (all interventional studies, or phase 2 or 3 only).
#'
#' @format A data.frame with columns `design`, `term`, `year`, `n`, `denominator`, `share_pct`, `download_date`.
#' @source ClinicalTrials.gov API v2, downloaded with `data-raw/ctgov_download.R`.
"ctgov_counts"

#' ClinicalTrials.gov interventional records mentioning "Bayesian"
#'
#' @format A data.frame with columns `nct`, `title`, `phase`, `sponsor_class`, `conditions`, `area` (keyword-mapped condition area), `start_date`, `start_year`, `status`, `enrollment`.
#' @source ClinicalTrials.gov API v2, downloaded with `data-raw/ctgov_download.R`.
"ctgov_bayesian_records"
