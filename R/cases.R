# Classical O'Brien-Fleming, 4 equally spaced looks, one-sided 0.025: gs_boundary(4)
.select_obf <- c(4.0487, 2.8628, 2.3375, 2.0243)
.select_fracs <- c(0.25, 0.5, 0.75, 1)

# SELECT-type reconstruction (Gaussian score process in information time).
# Convention: the effect is benefit-positive, theta = -log(HR), so theta > 0 means the treatment
# reduces the hazard. Information per event is 1/4, I = d/4, and the score S has mean theta * I;
# the z statistic is z = S / sqrt(I) and is positive for benefit.

#' Posterior probability of benefit
#'
#' Pr(theta > 0 | data) with benefit-positive `theta = -log(HR)`, prior N(0, tau^2), score S (mean `theta * I`, so S > 0 favours the treatment), information I.
#'
#' @param S score statistic (positive for benefit).
#' @param I information.
#' @param tau prior sd.
#' @return numeric.
#' @examples select_post_prob(10, 238)
#' @export
select_post_prob <- function(S, I, tau = 0.25) {
  A <- I + 1 / tau^2
  pnorm(S / sqrt(A))
}

#' Log Bayes factor (mixture e-value) for the score process
#'
#' Two-sided normal mixture, `theta ~ N(0, tau^2)` against `theta = 0`. It depends on `S^2`, so a rule using it for benefit must also require `S > 0`.
#'
#' @inheritParams select_post_prob
#' @return numeric.
#' @examples select_log_bf(10, 238)
#' @export
select_log_bf <- function(S, I, tau = 0.25) {
  A <- I + 1 / tau^2
  -0.5 * log(1 + I * tau^2) + S^2 / (2 * A)
}

#' Monte Carlo of the calibrated SELECT rules
#'
#' Simulates Gaussian score paths on a grid of event counts (steps of two, plus the interim, planned and total event counts), calibrates the posterior bar p* to type I error `alpha`, and compares four rules: a hypothetical O'Brien-Fleming design with four equally spaced looks, an e-value mixture, the Bayes rule, and SELECT's actual design. The effect is benefit-positive, `theta = -log(HR)`. Paths are processed in batches to bound memory.
#'
#' The actual SELECT design (Lincoff et al. 2023) is one superiority interim at two thirds of the planned events followed by the final analysis at the planned events, with O'Brien-Fleming-type boundaries at one-sided level `alpha`, `gs_boundary(2, alpha, t = c(interim_events / design_events, 1))`. The two other group sequential and the continuous rules run to `d_tot` events if they do not cross; the actual design ends at `design_events` if it does not cross.
#'
#' Both group sequential comparators use boundaries `gs_boundary(K, alpha, t = actual information fractions)`, so they respect `alpha` and the event rounding.
#'
#' Event summaries. For each rule, `median_events_success` and `mean_events_success` are conditional on crossing (among trials that cross the boundary at some look; these are discovery times, not trial sizes). `median_events_all` and `mean_events_all` are unconditional: trials that never cross are assigned the maximum event count of the rule (`d_tot`, or `design_events` for the actual design), so they are the median and expected number of events at which the trial ends. `power` is the probability of crossing under `HR`; `prob_cross` is the same number under an explicit name.
#'
#' @param d_tot total events (maximum for the four-look, e-value and Bayes rules).
#' @param HR true hazard ratio for the power calculation (benefit when below one).
#' @param nsim number of simulated trials (per hypothesis).
#' @param tau prior sd for the benefit-positive log hazard ratio `-log HR`.
#' @param alpha one-sided level.
#' @param seed optional seed.
#' @param batch trials per batch.
#' @param design_events planned events of SELECT's actual design (default `min(1225, d_tot)`).
#' @param interim_fraction information fraction of the actual design's interim (default 2/3).
#' @return list with `naive_typeI`, `pstar`, `bayes_typeI`, `ev_typeI`, `gs_typeI`, `actual_typeI`, `gs_bounds`, `actual_bounds`, `dg` (event grid), `max_events` (maximum event count per rule) and `rows`, a data.frame with `rule`, `median_events_success`, `mean_events_success`, `median_events_all`, `mean_events_all`, `power`, `prob_cross`.
#' @examples select_simulate(200, 0.8, nsim = 200, seed = 1)$pstar
#' @references Polson and Sokolov, case-study calibration for SELECT; Lincoff et al. (2023), NEJM 389, 2221-2232.
#' @export
select_simulate <- function(d_tot = 1270, HR = 0.80, nsim = 40000, tau = 0.25,
                            alpha = 0.025, seed = NULL, batch = 4000,
                            design_events = min(1225L, d_tot), interim_fraction = 2 / 3) {
  if (!is.null(seed)) set.seed(seed)
  if (design_events > d_tot) stop("design_events must not exceed d_tot")
  fracs <- .select_fracs
  theta <- -log(HR)
  d_int <- as.integer(round(interim_fraction * design_events))
  dg <- sort(unique(c(seq(4L, d_tot, by = 2L), as.integer(d_tot), d_int, as.integer(design_events))))
  Ig <- dg / 4
  K <- length(dg)
  dI <- diff(c(0, Ig))
  kk <- vapply(fracs, function(f) which.min(abs(dg - f * d_tot)), 1L)
  kk2 <- match(c(d_int, design_events), dg)
  B <- gs_boundary(length(kk), alpha, t = dg[kk] / dg[kk[length(kk)]])
  B2 <- gs_boundary(2, alpha, t = c(d_int / design_events, 1))
  sI <- sqrt(Ig)
  paths <- function(th, nb) {
    x <- matrix(rnorm(K * nb, sd = sqrt(dI)), K, nb)
    x <- apply(x, 2, cumsum)
    x + th * Ig
  }
  gs_hit <- function(S, k, bnd) {
    Z <- S[k, , drop = FALSE] / sI[k]
    colSums(Z >= bnd) > 0
  }
  starts <- seq(1, nsim, by = batch)
  sizes <- pmin(batch, nsim - starts + 1)
  # pass 1: null
  naive <- ev0 <- gs0 <- act0 <- logical(nsim)
  mx <- numeric(nsim)
  pos <- 0
  for (nb in sizes) {
    S <- paths(0, nb)
    pp <- select_post_prob(S, Ig, tau)
    idx <- pos + seq_len(nb)
    mx[idx] <- apply(pp, 2, max)
    ev0[idx] <- colSums(select_log_bf(S, Ig, tau) >= log(1 / alpha) & S > 0) > 0
    gs0[idx] <- gs_hit(S, kk, B)
    act0[idx] <- gs_hit(S, kk2, B2)
    pos <- pos + nb
  }
  kcut <- round(alpha * nsim)
  pstar <- sort(mx)[nsim - kcut + 1]
  # pass 2: alternative
  first <- list(gs = rep(NA_real_, nsim), ev = rep(NA_real_, nsim),
                bayes = rep(NA_real_, nsim), actual = rep(NA_real_, nsim))
  pos <- 0
  firstidx <- function(m) {
    hit <- colSums(m) > 0
    ifelse(hit, dg[apply(m, 2, which.max)], NA_real_)
  }
  for (nb in sizes) {
    S <- paths(theta, nb)
    idx <- pos + seq_len(nb)
    Z <- S / sI
    g <- matrix(FALSE, K, nb)
    for (j in seq_along(kk)) g[kk[j], ] <- Z[kk[j], ] >= B[j]
    first$gs[idx] <- firstidx(g)
    g2 <- matrix(FALSE, K, nb)
    for (j in seq_along(kk2)) g2[kk2[j], ] <- Z[kk2[j], ] >= B2[j]
    first$actual[idx] <- firstidx(g2)
    first$ev[idx] <- firstidx(select_log_bf(S, Ig, tau) >= log(1 / alpha) & S > 0)
    first$bayes[idx] <- firstidx(select_post_prob(S, Ig, tau) >= pstar)
    pos <- pos + nb
  }
  nm <- c(sprintf("O'Brien-Fleming, four scheduled looks, alpha = %g", alpha),
          sprintf("E-value mixture, continuous, E >= %g", 1 / alpha),
          sprintf("Bayes posterior, continuous, p* = %.3f", pstar),
          sprintf("SELECT actual design: interim at %d of %d events, then final", d_int, design_events))
  max_events <- c(d_tot, d_tot, d_tot, design_events)
  with_max <- function(x, m) ifelse(is.na(x), m, x)
  rows <- data.frame(
    rule = nm,
    median_events_success = vapply(first, median, 0, na.rm = TRUE),
    mean_events_success = vapply(first, mean, 0, na.rm = TRUE),
    median_events_all = unname(mapply(function(x, m) median(with_max(x, m)), first, max_events)),
    mean_events_all = unname(mapply(function(x, m) mean(with_max(x, m)), first, max_events)),
    power = vapply(first, function(x) mean(!is.na(x)), 0),
    row.names = NULL)
  rows$prob_cross <- rows$power
  # pass-1 naive error at the nominal 0.975 posterior bar
  naive_typeI <- mean(mx >= 0.975)
  list(naive_typeI = naive_typeI, pstar = pstar, bayes_typeI = mean(mx >= pstar),
       ev_typeI = mean(ev0), gs_typeI = mean(gs0), actual_typeI = mean(act0),
       gs_bounds = B, actual_bounds = B2, dg = dg, max_events = max_events, rows = rows)
}

#' Calendar month at which a fraction of events has accrued
#'
#' Linear enrolment over `T_en` months, events proportional to person time.
#'
#' @param d events at the decision.
#' @param d_tot total events.
#' @param T_en enrolment duration in months.
#' @param T_end trial duration in months.
#' @return month (numeric, vectorized over `d`).
#' @examples select_month_of_events(464, 1270)
#' @export
select_month_of_events <- function(d, d_tot, T_en = 36, T_end = 60) {
  pt <- function(t) ifelse(t <= T_en, 0.5 * t^2 / T_en, 0.5 * T_en + (t - T_en))
  one <- function(d) uniroot(function(t) pt(t) / pt(T_end) - d / d_tot,
                             c(1e-9, T_end), tol = 1e-10)$root
  vapply(d, one, 0)
}

#' Look table at the observed effect
#'
#' Benefit-positive convention: `theta = -log(HR)`, score `S = theta * I`, and the signed statistic `z = S / sqrt(I)`, which is positive for benefit (HR below one) and negative for harm.
#'
#' @param d_tot total events.
#' @param HR hazard ratio assumed observed at each look.
#' @param tau prior sd.
#' @return data.frame with `events`, `z` (signed, benefit-positive), `obf_bound`, `log_bf`, `post_prob` (probability of benefit).
#' @examples select_look_table(952, 0.8)
#' @export
select_look_table <- function(d_tot, HR, tau = 0.25) {
  f <- .select_fracs
  B <- .select_obf
  theta <- -log(HR)
  d <- f * d_tot
  I <- d / 4
  S <- theta * I
  data.frame(events = as.integer(round(d)), z = S / sqrt(I), obf_bound = B,
             log_bf = select_log_bf(S, I, tau), post_prob = select_post_prob(S, I, tau))
}
