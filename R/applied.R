# Applied examples: survival workbook, Gittins index, monitoring rulebooks, reanalyses.
# Effects are benefit-positive: theta = -log(HR) or -log(RR) unless stated otherwise.

# ---------------------------------------------------------------------------
# Survival
# ---------------------------------------------------------------------------

#' Events as information: posterior probabilities for a hazard ratio
#'
#' With one to one allocation and `D` events, `Var(log HR) = 4/D`, so the
#' information is `I = D/4`. Posterior probabilities for an observed hazard
#' ratio under a flat prior or a mean-zero normal prior on `theta = -log HR`.
#'
#' @param D Numbers of events.
#' @param hr Observed hazard ratio.
#' @param hr_mcid Hazard ratio defining a clinically meaningful reduction.
#' @param prior_sd Standard deviation of a N(0, prior_sd^2) prior on the log hazard ratio, or `NULL` for a flat prior.
#' @return data.frame with `D`, `se`, `p_benefit` (Pr(HR < 1)), `p_mcid` (Pr(HR < hr_mcid)) and `prior_events` (`4/prior_sd^2`, 0 if flat).
#' @references Schoenfeld (1981) Biometrika 68, 316-319.
#' @examples
#' surv_events_table()
#' @export
surv_events_table <- function(D = c(100, 200, 400, 800), hr = 0.8, hr_mcid = 0.9, prior_sd = NULL) {
  se <- 2 / sqrt(D)
  theta <- -log(hr)
  prec <- 1 / se^2 + if (is.null(prior_sd)) 0 else 1 / prior_sd^2
  s1 <- sqrt(1 / prec)
  m1 <- (theta / se^2) / prec
  data.frame(D = D, se = se,
             p_benefit = stats::pnorm(m1 / s1),
             p_mcid = stats::pnorm((m1 + log(hr_mcid)) / s1),
             prior_events = if (is.null(prior_sd)) 0 else 4 / prior_sd^2)
}

#' Conditional power and predictive probability in an event-driven trial
#'
#' A trial with `Dmax` events succeeds if the final one-sided z exceeds `z_crit`.
#' At an interim with `D1` events and observed hazard ratio `hr_obs`, returns the
#' interim z, the flat-prior predictive probability of success, and conditional
#' power under the design hazard ratio and under the current trend.
#'
#' @param D1 Interim events (vector).
#' @param hr_obs Observed hazard ratio at the interim (vector, recycled).
#' @param Dmax Planned total events.
#' @param hr_design Design hazard ratio.
#' @param z_crit Final success threshold on the benefit-positive z scale.
#' @return data.frame with `D1`, `hr_obs`, `z`, `pp`, `cp_design`, `cp_trend`.
#' @examples
#' surv_pp_cp(c(300, 300, 300, 450), c(0.75, 0.85, 0.95, 0.88))
#' @export
surv_pp_cp <- function(D1, hr_obs, Dmax = 600, hr_design = 0.8, z_crit = 1.96) {
  n <- max(length(D1), length(hr_obs))
  D1 <- rep_len(D1, n); hr_obs <- rep_len(hr_obs, n)
  D2 <- Dmax - D1
  th <- -log(hr_obs)
  crit <- z_crit * 2 / sqrt(Dmax)           # final estimate must exceed this
  sd_cp <- (D2 / Dmax) * 2 / sqrt(D2)
  v_pp <- (D2 / Dmax)^2 * (4 / D2 + 4 / D1)
  data.frame(D1 = D1, hr_obs = hr_obs, z = th / (2 / sqrt(D1)),
             pp = stats::pnorm((th - crit) / sqrt(v_pp)),
             cp_design = stats::pnorm(((D1 * th + D2 * (-log(hr_design))) / Dmax - crit) / sd_cp),
             cp_trend = stats::pnorm((th - crit) / sd_cp))
}

#' Gamma posterior for a constant hazard
#'
#' Posterior Gamma(a + d, b + T) for `d` events in `T` patient-years, with the
#' posterior mean and equal-tailed interval for the hazard and for the
#' exponential median survival `log 2 / lambda` (exact, by the reciprocal transformation).
#'
#' @param d Events. @param T Patient-years.
#' @param a,b Gamma prior shape and rate (vague default).
#' @param level Interval level.
#' @return Named list: `hazard` (mean, lo, hi), `median` (mean by Monte Carlo-free closed form, lo, hi) and `shape`, `rate`.
#' @examples
#' surv_gamma_hazard(60, 400)
#' @export
surv_gamma_hazard <- function(d, T, a = 0.01, b = 0.01, level = 0.95) {
  sh <- a + d; rt <- b + T
  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  haz <- stats::qgamma(q, sh, rt)
  # median = log2 / lambda ; E[1/lambda] = rt / (sh - 1)
  list(hazard = c(mean = sh / rt, lo = haz[1], hi = haz[2]),
       median = c(mean = log(2) * rt / (sh - 1), lo = log(2) / haz[2], hi = log(2) / haz[1]),
       shape = sh, rate = rt)
}

#' Survival function of a piecewise-exponential distribution
#'
#' @param t Times. @param rates Hazards in each interval. @param breaks Left endpoints of the intervals (first is 0).
#' @return Survival probabilities `S(t)`.
#' @examples
#' surv_pwexp_S(c(1, 5), c(0.15, 0.1), c(0, 2))
#' @export
surv_pwexp_S <- function(t, rates, breaks = 0) {
  stopifnot(length(rates) == length(breaks), breaks[1] == 0)
  ends <- c(breaks[-1], Inf)
  vapply(t, function(tt) exp(-sum(rates * pmax(0, pmin(tt, ends) - breaks))), numeric(1))
}

#' Restricted mean survival time for a piecewise-exponential arm
#'
#' Closed form: on each interval, `S` decays exponentially from its value at the
#' left endpoint, so the integral is `S(b_k)(1 - exp(-lambda_k w_k))/lambda_k`.
#' With a cure fraction `pi`, `S(t) = pi + (1 - pi) S0(t)`.
#'
#' @param tau Horizon. @inheritParams surv_pwexp_S
#' @param cure Cure fraction (default 0).
#' @return RMST(tau).
#' @examples
#' surv_rmst(5, 0.15)
#' @export
surv_rmst <- function(tau, rates, breaks = 0, cure = 0) {
  ends <- pmin(c(breaks[-1], Inf), tau)
  w <- pmax(0, ends - breaks)
  Sb <- vapply(breaks, function(b) surv_pwexp_S(b, rates, breaks), numeric(1))
  r0 <- sum(Sb * ifelse(rates > 0, (1 - exp(-rates * w)) / rates, w))
  cure * tau + (1 - cure) * r0
}

#' Non-proportional hazards examples: delayed effect and cure model
#'
#' Control arm exponential with hazard `h0`. Delayed effect: treated hazard
#' equals `h0` before `delay` and `hr * h0` afterwards. Cure model: a fraction
#' `cure` of treated patients never fail, the rest follow the control hazard.
#' Returns RMST by arm and its difference at the horizon, and the
#' average hazard ratio `log S1(tau)/log S0(tau)` (the ratio of cumulative
#' hazards, constant only under proportional hazards).
#'
#' @param tau Horizon (years). @param h0 Control hazard. @param delay Delay in years.
#' @param hr Hazard ratio after the delay. @param cure Cure fraction in the cure-model scenario.
#' @return data.frame with one row per scenario: `scenario`, `rmst_control`, `rmst_treated`, `rmst_diff`, `cumhaz_ratio`.
#' @examples
#' surv_nph_example()
#' @export
surv_nph_example <- function(tau = 5, h0 = 0.15, delay = 1, hr = 0.6, cure = 0.2) {
  r0 <- surv_rmst(tau, h0)
  row <- function(lab, r1, S1) data.frame(scenario = lab, rmst_control = r0, rmst_treated = r1,
                                          rmst_diff = r1 - r0, cumhaz_ratio = log(S1) / log(exp(-h0 * tau)))
  rbind(
    row("proportional", surv_rmst(tau, h0 * hr), exp(-h0 * hr * tau)),
    row("delayed effect", surv_rmst(tau, c(h0, h0 * hr), c(0, delay)),
        surv_pwexp_S(tau, c(h0, h0 * hr), c(0, delay))),
    row("cure model", surv_rmst(tau, h0, cure = cure), cure + (1 - cure) * exp(-h0 * tau)))
}

# ---------------------------------------------------------------------------
# Gittins index
# ---------------------------------------------------------------------------

#' Gittins index of a Bernoulli arm with a Beta(a, b) state
#'
#' Calibration against a standard arm paying `lambda` per period forever
#' (retirement reward `lambda/(1 - gamma)`), solving the one-armed problem by
#' backward induction on the states reachable from Beta(a, b) within `N` further
#' pulls (terminal value `max(R, p/(1-gamma))`), and bisecting on `lambda`.
#' Truncation error is of order `gamma^N`.
#'
#' @param a,b Beta parameters of the current posterior.
#' @param gamma Discount factor in (0, 1).
#' @param N Depth of the backward induction.
#' @param tol Bisection tolerance.
#' @return The Gittins index (a number in (0, 1)).
#' @references Gittins and Jones (1974); Berry and Fristedt (1985).
#' @examples
#' gittins_index(1, 1, 0.9, N = 100)
#' @export
gittins_index <- function(a, b, gamma = 0.9, N = 150, tol = 1e-6) {
  stopifnot(gamma > 0, gamma < 1, a > 0, b > 0)
  index_for <- function(lam) {
    R <- lam / (1 - gamma)
    s <- 0:N
    V <- pmax(R, (a + s) / (a + b + N) / (1 - gamma))
    for (n in (N - 1):0) {
      s <- 0:n; p <- (a + s) / (a + b + n)
      cont <- p * (1 + gamma * V[s + 2]) + (1 - p) * gamma * V[s + 1]
      V <- pmax(R, cont)
    }
    V[1] > R + 1e-12
  }
  lo <- 0; hi <- 1
  while (hi - lo > tol) {
    m <- (lo + hi) / 2
    if (index_for(m)) lo <- m else hi <- m
  }
  (lo + hi) / 2
}

#' Gittins index table
#'
#' @param states data.frame with columns `a`, `b` (default: the five states of the rl chapter).
#' @param gamma Vector of discount factors.
#' @param N Depth; if `NULL`, `ceiling(log(1e-4)/log(gamma))` capped at 500.
#' @return data.frame with `a`, `b`, `mean`, and one column `gamma_<g>` per discount factor.
#' @examples
#' gittins_table(gamma = 0.9, N = 100)
#' @export
gittins_table <- function(states = data.frame(a = c(1, 2, 1, 5, 10), b = c(1, 1, 2, 5, 10)),
                          gamma = c(0.9, 0.99), N = NULL) {
  out <- data.frame(a = states$a, b = states$b, mean = states$a / (states$a + states$b))
  for (g in gamma) {
    n <- if (is.null(N)) min(500, ceiling(log(1e-4) / log(g))) else N
    out[[paste0("gamma_", g)]] <- mapply(gittins_index, states$a, states$b, MoreArgs = list(gamma = g, N = n))
  }
  out
}

# ---------------------------------------------------------------------------
# Monitoring: three rulebooks plus a calibrated Bayes rule
# ---------------------------------------------------------------------------

# z-boundary (benefit-positive) at which the posterior Pr(theta > 0) under a N(0, s0^2) prior reaches `thr`
.monitor_post_bound <- function(thr, I, prior_sd) {
  # posterior prob = Phi(S / sqrt(I + 1/s0^2)) with S = z sqrt(I)
  stats::qnorm(thr) * sqrt(I + 1 / prior_sd^2) / sqrt(I)
}

#' Calibrate a posterior-probability threshold to a one-sided error rate
#'
#' The Bayes rule "stop for efficacy at the first look where
#' Pr(theta > 0 | data) >= thr" under a N(0, prior_sd^2) prior (flat if `prior_sd = Inf`)
#' is a z-boundary `z_k = qnorm(thr) sqrt(1 + 1/(I_k prior_sd^2))`. The null
#' maximum over the looks is computed exactly with [gs_crossing_prob()] and `thr`
#' solved so that it equals `alpha`.
#'
#' @param events Events at the looks. @param prior_sd Prior sd on the log hazard ratio (`Inf` = flat).
#' @param alpha One-sided type I error.
#' @return list with `threshold`, `bounds` (z boundaries) and `alpha`.
#' @examples
#' monitor_calibrate_bayes(c(300, 600, 900, 1200), prior_sd = 0.1)
#' @export
monitor_calibrate_bayes <- function(events, prior_sd = Inf, alpha = 0.025) {
  I <- events / 4; t <- events / max(events)
  bnd <- function(thr) .monitor_post_bound(thr, I, prior_sd)
  f <- function(x) gs_crossing_prob(bnd(stats::pnorm(x)), t)$total - alpha
  x <- stats::uniroot(f, c(1, 8), tol = 1e-10)$root
  list(threshold = stats::pnorm(x), bounds = bnd(stats::pnorm(x)), alpha = alpha)
}

#' Three rulebooks, side by side, for an outcomes trial (plus calibrated Bayes)
#'
#' For observed hazard ratios at each look reports the z statistic, the
#' O'Brien-Fleming boundary, flat and skeptical posterior probabilities of
#' benefit and of a reduction beyond `hr_mcid`, the anytime-valid Bayes factor
#' (N(0, e_sd^2) mixture against the null, evaluated as the normal marginal
#' likelihood ratio at the observed estimate, as in the zip), and the
#' calibrated Bayes rule: the skeptical posterior probability against a threshold
#' chosen so the null probability of ever exceeding it over the looks is `alpha`.
#'
#' @param events Events at the looks. @param hr Observed hazard ratios.
#' @param skeptic_sd Skeptical prior sd. @param e_sd Mixture sd for the Bayes factor.
#' @param hr_mcid Hazard ratio of a meaningful reduction.
#' @param obf_exact If `TRUE` use [gs_boundary()] (OBF constant solved exactly) else the same; kept for symmetry.
#' @param alpha One-sided type I error.
#' @return data.frame with `events`, `hr`, `z`, `obf_bound`, `obf_stop`, `flat_p`, `flat_p_mcid`, `skep_p`, `skep_p_mcid`,
#'   `bf`, `bf_stop` (bf >= 1/alpha), `cal_threshold`, `cal_bound`, `cal_stop`. The calibrated threshold is also in `attr(, "calibration")`.
#' @examples
#' monitor_three_rulebooks()
#' @export
monitor_three_rulebooks <- function(events = c(300, 600, 900, 1200), hr = c(0.78, 0.83, 0.84, 0.85),
                                    skeptic_sd = 0.1, e_sd = 0.2, hr_mcid = 0.9,
                                    obf_exact = TRUE, alpha = 0.025) {
  I <- events / 4; se <- 2 / sqrt(events); t <- events / max(events)
  th <- -log(hr); z <- th / se
  obf <- gs_boundary(length(events), alpha, "obf", t = t)
  post <- function(sd0) {
    prec <- 1 / se^2 + if (is.finite(sd0)) 1 / sd0^2 else 0
    m <- (th / se^2) / prec; s <- sqrt(1 / prec)
    list(p = stats::pnorm(m / s), pm = stats::pnorm((m + log(hr_mcid)) / s), m = m)
  }
  fl <- post(Inf); sk <- post(skeptic_sd)
  bf <- stats::dnorm(th, 0, sqrt(se^2 + e_sd^2)) / stats::dnorm(th, 0, se)
  cal <- monitor_calibrate_bayes(events, skeptic_sd, alpha)
  out <- data.frame(events = events, hr = hr, z = z, obf_bound = obf, obf_stop = z >= obf,
                    flat_p = fl$p, flat_p_mcid = fl$pm, skep_p = sk$p, skep_p_mcid = sk$pm,
                    skep_mean_hr = exp(-sk$m), bf = bf, bf_stop = bf >= 1 / alpha,
                    cal_threshold = cal$threshold, cal_bound = cal$bounds,
                    cal_stop = sk$p >= cal$threshold)
  attr(out, "calibration") <- cal
  out
}

#' Harm and futility summaries for an outcomes-trial monitor
#'
#' Flat-prior posterior probability that the hazard ratio exceeds `hr_harm`
#' (Pr(HR > hr_harm)), the observed hazard ratio at which that probability
#' equals `p_stop` (which is `hr_harm` itself when `p_stop = 0.5`), and the
#' predictive probability of final success (futility summary, [surv_pp_cp()]).
#'
#' @param events Events at the look. @param hr Observed hazard ratio.
#' @param hr_harm Harm hazard ratio. @param p_stop Harm stopping threshold.
#' @param Dmax Planned events. @param z_crit Final z threshold.
#' @return Named list `p_harm`, `trigger_hr`, `pp_success`.
#' @examples
#' monitor_harm_rule(60, 1.3)
#' monitor_harm_rule(600, 0.83, Dmax = 1200, z_crit = 2.02)
#' @export
monitor_harm_rule <- function(events, hr, hr_harm = 1.2, p_stop = 0.5, Dmax = events, z_crit = 1.96) {
  se <- 2 / sqrt(events)
  p_harm <- stats::pnorm((log(hr) - log(hr_harm)) / se)
  trigger <- exp(log(hr_harm) + stats::qnorm(p_stop) * se)
  pp <- if (events < Dmax) surv_pp_cp(events, hr, Dmax = Dmax, z_crit = z_crit)$pp else NA_real_
  list(p_harm = p_harm, trigger_hr = trigger, pp_success = pp)
}

# ---------------------------------------------------------------------------
# Reanalyses
# ---------------------------------------------------------------------------

#' Normal-approximation posterior with thresholds
#'
#' Posterior for a benefit-positive effect `theta` from estimate `est`,
#' standard error `se` and a normal prior (flat if `prior_sd = Inf`).
#'
#' @param est,se Estimate and standard error (benefit-positive scale).
#' @param prior_mean,prior_sd Prior mean and sd (`Inf` for flat).
#' @param mcid Benefit threshold(s) on the theta scale.
#' @return data.frame with `post_mean`, `post_sd`, `p_benefit` (Pr(theta > 0)) and one column `p_gt_<mcid>` per threshold.
#' @examples
#' reanalysis_normal(0.01, 0.003, 0, 0.0035, mcid = 0.01)
#' @export
reanalysis_normal <- function(est, se, prior_mean = 0, prior_sd = Inf, mcid = numeric(0)) {
  pp <- if (is.finite(prior_sd)) 1 / prior_sd^2 else 0
  prec <- 1 / se^2 + pp
  s1 <- sqrt(1 / prec)
  m1 <- (est / se^2 + if (pp > 0) prior_mean * pp else 0) / prec
  out <- data.frame(post_mean = m1, post_sd = s1, p_benefit = stats::pnorm(m1 / s1))
  for (m in mcid) out[[paste0("p_gt_", format(m))]] <- stats::pnorm((m1 - m) / s1)
  out
}

#' Reanalysis of a ratio from its published point estimate and 95% interval
#'
#' Works on `theta = -log(ratio)` (benefit positive for HR, RR, OR below one),
#' `se = (log u - log l)/(2 * qnorm(0.975))`, and several priors given on the
#' log-ratio scale as `list(name = c(mean, sd))` with `sd = Inf` for flat.
#'
#' @param ratio,lo,hi Published ratio and 95% interval. @param se Optional override of the standard error.
#' @param priors Named list of `c(mean_of_log_ratio, sd)`.
#' @param thresholds Ratios (e.g. `c(0.8, 0.67)`) for Pr(ratio < r).
#' @return data.frame with `prior`, `post_ratio`, `p_benefit` (Pr(ratio < 1)) and `p_<r>` columns; `se` and `z` in attributes.
#' @examples
#' reanalysis_ratio(0.75, 0.55, 1.02, list(flat = c(0, Inf), skeptic = c(0, 0.15)), 0.8)
#' @export
reanalysis_ratio <- function(ratio, lo, hi, priors = list(flat = c(0, Inf)), thresholds = numeric(0), se = NULL) {
  if (is.null(se)) se <- (log(hi) - log(lo)) / (2 * stats::qnorm(0.975))
  rows <- lapply(names(priors), function(nm) {
    pr <- priors[[nm]]
    r <- reanalysis_normal(-log(ratio), se, -pr[1], pr[2], mcid = -log(thresholds))
    out <- data.frame(prior = nm, post_ratio = exp(-r$post_mean), p_benefit = r$p_benefit)
    for (k in seq_along(thresholds)) out[[paste0("p_", format(thresholds[k]))]] <- r[[3 + k]]
    out
  })
  out <- do.call(rbind, rows)
  attr(out, "se") <- se
  out
}

#' GUSTO prior-weight table (Brophy and Joseph style)
#'
#' Risk difference (streptokinase minus tPA, benefit positive for tPA) from the
#' reported arm sizes and event rates, combined with a prior centred at zero,
#' with sd `prior_se / sqrt(w)`, representing GISSI-2 and ISIS-3 at weight `w`.
#'
#' @param n_tpa,n_sk Arm sizes. @param p_tpa,p_sk Mortality.
#' @param prior_se Prior standard error at full weight. @param w Prior weights (`0` = flat).
#' @param mcid Clinically important risk difference.
#' @return data.frame with `prior`, `post_mean`, `p_benefit`, `p_mcid`; `attr(, "est")` and `attr(, "se")`.
#' @examples
#' reanalysis_gusto()
#' @export
reanalysis_gusto <- function(n_tpa = 10400, n_sk = 20200, p_tpa = 0.063, p_sk = 0.073,
                             prior_se = 0.0035, w = c(0, 0.1, 0.5, 1), mcid = 0.01) {
  est <- p_sk - p_tpa
  se <- sqrt(p_tpa * (1 - p_tpa) / n_tpa + p_sk * (1 - p_sk) / n_sk)
  rows <- lapply(w, function(wi) {
    r <- reanalysis_normal(est, se, 0, if (wi == 0) Inf else prior_se / sqrt(wi), mcid = mcid)
    data.frame(prior = if (wi == 0) "Flat" else paste0("w = ", wi), post_mean = r$post_mean,
               p_benefit = r$p_benefit, p_mcid = r[[4]])
  })
  out <- do.call(rbind, rows)
  attr(out, "est") <- est; attr(out, "se") <- se
  out
}

#' EOLIA reanalysis (relative risk of death at 60 days, from counts)
#'
#' @param y1,n1 Deaths and size, ECMO. @param y0,n0 Deaths and size, control.
#' @return data.frame as [reanalysis_ratio()] with thresholds 0.8 and 0.67; `attr(, "se")`.
#' @examples
#' reanalysis_eolia()
#' @export
reanalysis_eolia <- function(y1 = 44, n1 = 124, y0 = 57, n0 = 125) {
  lrr <- log((y1 / n1) / (y0 / n0)); se <- sqrt(1 / y1 - 1 / n1 + 1 / y0 - 1 / n0)
  reanalysis_ratio(exp(lrr), 0, 0, se = se,
                   priors = list(Flat = c(0, Inf), "Skeptical, s0 = 0.25" = c(0, 0.25),
                                 "Skeptical, s0 = 0.15" = c(0, 0.15),
                                 Enthusiastic = c(log(0.67), 0.25)),
                   thresholds = c(0.8, 0.67))
}

#' ANDROMEDA SHOCK reanalysis (hazard ratio 0.75, 95% CI 0.55 to 1.02)
#'
#' @return data.frame as [reanalysis_ratio()] with threshold 0.8.
#' @examples
#' reanalysis_andromeda()
#' @export
reanalysis_andromeda <- function() {
  reanalysis_ratio(0.75, 0.55, 1.02,
                   priors = list(Flat = c(0, Inf), "Neutral, s0 = 0.35" = c(0, 0.35),
                                 "Skeptical, s0 = 0.15" = c(0, 0.15),
                                 "Pessimistic, centered at HR = 1.1" = c(log(1.1), 0.15),
                                 "Optimistic, centered at HR = 0.8" = c(log(0.8), 0.2)),
                   thresholds = 0.8)
}

#' Hierarchical shrinkage of stratum effects (uniform prior on tau, flat on mu)
#'
#' Normal hierarchical model `est_j ~ N(theta_j, se_j^2)`, `theta_j ~ N(mu, tau^2)`,
#' flat prior on `mu`, uniform prior on `tau` over `(0, tau_max]` (grid). Posterior
#' means and probabilities of benefit are averaged over the posterior of `tau`
#' (with `mu` integrated out, so the conditional variance includes the variance of `mu`).
#' Estimates are on the benefit-positive scale.
#'
#' @param est,se Stratum estimates and standard errors. @param tau_max Upper bound for tau. @param ngrid Grid size.
#' @return list: `strata` data.frame (`est`, `se`, `shrunk`, `p_benefit`, `p_benefit_flat`), `tau_mean`, `tau_grid`, `tau_post`.
#' @examples
#' reanalysis_hier(c(0.44, 0.2, -0.17), c(0.12, 0.07, 0.13))$tau_mean
#' @export
reanalysis_hier <- function(est, se, tau_max = 1, ngrid = 400) {
  taus <- seq(0.001, tau_max, length.out = ngrid)
  ll <- vapply(taus, function(t) {
    v <- se^2 + t^2; w <- 1 / v; mu <- sum(w * est) / sum(w)
    -0.5 * sum(log(v) + (est - mu)^2 / v) + 0.5 * log(1 / sum(w))
  }, numeric(1))
  pw <- exp(ll - max(ll)); pw <- pw / sum(pw)
  m <- rep(0, length(est)); pb <- m
  for (i in seq_along(taus)) {
    t <- taus[i]; v <- se^2 + t^2; w <- 1 / v; mu <- sum(w * est) / sum(w); vm <- 1 / sum(w)
    B <- se^2 / v; mj <- B * mu + (1 - B) * est; sj <- sqrt(se^2 * (1 - B) + B^2 * vm)
    m <- m + pw[i] * mj; pb <- pb + pw[i] * stats::pnorm(mj / sj)
  }
  list(strata = data.frame(est = est, se = se, shrunk = m, p_benefit = pb,
                           p_benefit_flat = stats::pnorm(est / se)),
       tau_mean = sum(pw * taus), tau_grid = taus, tau_post = pw)
}

#' RECOVERY dexamethasone subgroup shrinkage by respiratory support
#'
#' Published rate ratios and 95% intervals for invasive ventilation, oxygen only
#' and no oxygen, shrunk with [reanalysis_hier()].
#'
#' @param tau_max Upper bound of the uniform prior on tau.
#' @return data.frame with `stratum`, `rr`, `shrunk_rr`, `p_benefit` (Pr(RR < 1)), `p_benefit_flat`; `attr(, "tau_mean")`.
#' @examples
#' reanalysis_recovery()
#' @export
reanalysis_recovery <- function(tau_max = 1) {
  rr <- c(0.64, 0.82, 1.19); lo <- c(0.51, 0.72, 0.92); hi <- c(0.81, 0.94, 1.55)
  se <- (log(hi) - log(lo)) / (2 * stats::qnorm(0.975))
  h <- reanalysis_hier(-log(rr), se, tau_max)
  out <- data.frame(stratum = c("Invasive ventilation", "Oxygen only", "No oxygen"), rr = rr,
                    shrunk_rr = exp(-h$strata$shrunk), p_benefit = h$strata$p_benefit,
                    p_benefit_flat = h$strata$p_benefit_flat)
  attr(out, "tau_mean") <- h$tau_mean
  out
}
