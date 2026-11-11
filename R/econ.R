#' Modeled exposure after an early decision
#'
#' Converts the expected number of events at an early decision into the
#' exposure, placebo events and excess placebo events that would come after it.
#' This is a modeled counterfactual: events are assumed uniform in person time,
#' placebo participants are assumed switched to active treatment at the decision,
#' and the decision is assumed acted on.
#'
#' @param mean_events Expected events at the decision, averaged over all simulated trials (a trial that never crosses counts at its maximum).
#' @param d_tot Total events at the end of the trial.
#' @param n Participants randomized (1:1).
#' @param years Mean follow-up in years.
#' @param ev_placebo,ev_active Observed primary events in the placebo and active arms.
#' @return A list with `share` (expected share of events at the decision), `placebo_py`, `placebo_events` and `excess` (the modeled number of placebo-arm primary events in excess of the active arm's rate that switching placebo participants to active treatment at the decision would avert; it is not the number of placebo events remaining, which is `placebo_events`).
#' @examples exposure_after_decision(900, 1270, 17604, 39.8 / 12, 701, 569)
#' @export
exposure_after_decision <- function(mean_events, d_tot, n, years, ev_placebo, ev_active) {
  share <- mean_events / d_tot
  list(share = share, placebo_py = (1 - share) * n * years / 2,
       placebo_events = (1 - share) * ev_placebo, excess = (1 - share) * (ev_placebo - ev_active))
}

#' Modeled trial cost saving from an early decision
#'
#' The saving is the avoided share of exposure times the variable share of cost times the total cost. The variable share and the total cost are assumptions.
#'
#' @param share Expected share of events (and person time) accrued at the decision.
#' @param total Total trial cost.
#' @param v Variable share of cost, between 0 and 1.
#' @return The saving, in the units of `total`.
#' @examples trial_cost_saving(0.7, 157, 0.6)
#' @export
trial_cost_saving <- function(share, total, v) (1 - share) * v * total

#' Label-driven revenue brought forward by an earlier indication
#'
#' With patent expiry fixed, an indication granted `years` earlier adds `years` of label-driven sales: annual sales times the share attributable to the label times the years. Revenue, not profit, and an upper bound, since approval, reimbursement and uptake would not all move forward by the full statistical acceleration. The label share is an assumption.
#'
#' @param sales Annual sales of the drug (any currency unit).
#' @param label_share Share of sales attributable to the indication (scalar or vector).
#' @param years Years the indication is brought forward.
#' @return Revenue brought forward, in the units of `sales`.
#' @examples label_revenue_forward(8.4, c(0.10, 0.25), 2.4)
#' @export
label_revenue_forward <- function(sales, label_share, years) sales * label_share * years
