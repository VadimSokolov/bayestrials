# Safety, phase II monitoring, BLRM and ECMO tables (zip chapters safety, phase2, dose, ecmo).

# ---------------------------------------------------------------- rofecoxib

#' Cumulative rofecoxib posterior table
#'
#' Normal approximation on the log relative risk scale: the likelihood is
#' \eqn{N(\hat\theta, s^2)} and the prior is \eqn{N(0, s_0^2)} (or flat).
#' Reproduces the rofecoxib table of the safety chapter.
#'
#' @param rr,lo,hi Pooled relative risk and 95 percent interval limits, one per period.
#' @param period Labels for the periods.
#' @param s0 Prior standard deviations on the log RR scale (`Inf` = flat).
#' @param thresholds Relative-risk thresholds for the exceedance probabilities.
#' @return Data frame with period, prior, posterior RR (exp of the posterior mean
#'   of log RR), and `Pr(RR > t)` columns.
#' @references Juni P et al. (2004). Lancet 364:2021-2029.
#' @examples
#' safety_rofecoxib()
#' @export
safety_rofecoxib <- function(rr = c(2.30, 2.24), lo = c(1.22, 1.24), hi = c(4.33, 4.02),
                             period = c("End 2000", "End 2001"),
                             s0 = c(Inf, 0.3, 0.2), thresholds = c(1, 1.5)) {
  out <- list()
  for (i in seq_along(rr)) {
    th <- log(rr[i]); se <- (log(hi[i]) - log(lo[i])) / (2 * stats::qnorm(0.975))
    if (i == 1) { th <- 0.833; se <- 0.323 }  # zip rounds the end-2000 inputs
    for (p in s0) {
      v <- if (is.finite(p)) 1 / (1 / se^2 + 1 / p^2) else se^2
      m <- if (is.finite(p)) v * th / se^2 else th
      pr <- vapply(thresholds, function(t) 1 - stats::pnorm((log(t) - m) / sqrt(v)), 0)
      d <- data.frame(period = period[i], prior = if (is.finite(p)) paste0("skeptical s0=", p) else "flat",
                      post_rr = exp(m), log_rr = th, se = se)
      for (j in seq_along(thresholds)) d[[paste0("pr_rr_gt_", thresholds[j])]] <- pr[j]
      out[[length(out) + 1]] <- d
    }
  }
  do.call(rbind, out)
}

#' Bayes factor and posterior-probability monitoring for a safety signal
#'
#' Bayes factor of \eqn{N(0, \tau^2)} alternative against the null for a normal
#' estimate (the rofecoxib number 5.7), the Ville threshold \eqn{1/\alpha}, and,
#' optionally, a simulated false-alarm comparison under the null over repeated
#' looks with a posterior-probability rule and a Bayes-factor rule.
#'
#' @param est,se Estimate and standard error on the log RR scale.
#' @param tau Prior sd under the alternative.
#' @param alpha Ville level; the BF threshold is `1/alpha`.
#' @param looks Number of monthly looks in the simulation (0 = skip).
#' @param pp_threshold Posterior-probability threshold, Pr(theta > 0) with a N(0, tau^2) prior.
#' @param n_per_look Information accrued per look (precision of theta is `n_per_look * t`).
#' @param nsim Simulated null trials.
#' @param seed Random seed.
#' @return Named list: `bf`, `threshold`, `flagged`, and `false_alarm` (named vector, or NULL).
#' @examples
#' safety_bf_monitor(0.833, 0.323, 0.5)$bf
#' @export
safety_bf_monitor <- function(est = 0.833, se = 0.323, tau = 0.5, alpha = 0.05,
                              looks = 0, pp_threshold = 0.95, n_per_look = 4,
                              nsim = 20000, seed = 1) {
  bf <- stats::dnorm(est, 0, sqrt(se^2 + tau^2)) / stats::dnorm(est, 0, se)
  fa <- NULL
  if (looks > 0) {
    set.seed(seed)
    # Brownian information scale: at look t, precision = n_per_look * t, S_t = sum of N(0, 1/n_per_look)
    S <- t(apply(matrix(stats::rnorm(nsim * looks), nsim), 1, cumsum)) / sqrt(n_per_look)
    info <- n_per_look * seq_len(looks)
    # theta_hat_t = S_t * sqrt(n_per_look) / info  with S_t built from unit-variance steps scaled so Var(theta_hat_t)=1/info
    that <- sweep(S * n_per_look, 2, info, "/")
    prec <- sweep(matrix(0, nsim, looks), 2, info, "+")
    v <- 1 / (prec + 1 / tau^2); m <- v * prec * that
    pp <- 1 - stats::pnorm(-m / sqrt(v))
    bfs <- sqrt(1 / (1 + tau^2 * prec)) * exp(0.5 * prec^2 * that^2 / (prec + 1 / tau^2))
    fa <- c(posterior_prob = mean(apply(pp >= pp_threshold, 1, any)),
            bayes_factor = mean(apply(bfs >= 1 / alpha, 1, any)),
            ville_bound = alpha)
  }
  list(bf = bf, threshold = 1 / alpha, flagged = bf >= 1 / alpha, false_alarm = fa)
}

# ---------------------------------------------------------------- gamma-Poisson shrinker

.gps_logmarg <- function(n, E, a, b) {
  # negative binomial: size a, prob b/(b+E)
  stats::dnbinom(n, size = a, prob = b / (b + E), log = TRUE)
}

#' Gamma-Poisson shrinker (DuMouchel 1999) posterior summaries
#'
#' For counts `n` with expected counts `E`, a two-component gamma mixture prior on
#' the relative reporting ratio gives a two-component gamma mixture posterior,
#' component k being Gamma(alpha_k + n, beta_k + E).
#'
#' @param n,E Observed and expected counts.
#' @param P Prior weight of component 1 (the "signal" component).
#' @param a1,b1,a2,b2 Gamma shape and rate of the two components.
#' @param prob Quantile level for the lower bound (EB05 = 0.05).
#' @return Data frame: n, E, ratio, weight (posterior weight on component 1),
#'   ebgm (geometric mean exp(E log lambda)), eb05 (the `prob` quantile).
#' @references DuMouchel W (1999). Bayesian data mining in large frequency tables. The American Statistician 53:177-190.
#' @examples
#' safety_gps(c(1, 3, 12, 40, 2), c(0.05, 0.5, 3, 10, 2))
#' @export
safety_gps <- function(n, E, P = 1/3, a1 = 0.2, b1 = 0.1, a2 = 2, b2 = 4, prob = 0.05) {
  stopifnot(length(n) == length(E))
  l1 <- log(P) + .gps_logmarg(n, E, a1, b1)
  l2 <- log(1 - P) + .gps_logmarg(n, E, a2, b2)
  w <- 1 / (1 + exp(l2 - l1))
  A1 <- a1 + n; B1 <- b1 + E; A2 <- a2 + n; B2 <- b2 + E
  elog <- w * (digamma(A1) - log(B1)) + (1 - w) * (digamma(A2) - log(B2))
  q <- vapply(seq_along(n), function(i) {
    f <- function(x) w[i] * stats::pgamma(x, A1[i], B1[i]) + (1 - w[i]) * stats::pgamma(x, A2[i], B2[i]) - prob
    stats::uniroot(f, c(0, 1e3), tol = 1e-10)$root
  }, 0)
  data.frame(n = n, E = E, ratio = n / E, weight = w, ebgm = exp(elog), eb05 = q)
}

#' Fit the gamma-Poisson shrinker hyperparameters
#'
#' Maximum marginal likelihood for the two-component gamma mixture, by `optim`
#' (L-BFGS-B on log scale parameters and logit weight) from several starts.
#'
#' @param n,E Count and expected-count vectors.
#' @param start Starting values `c(P, a1, b1, a2, b2)`.
#' @return List: `par` (named P, a1, b1, a2, b2), `loglik`, `convergence`.
#' @examples
#' set.seed(1); E <- runif(300, 0.2, 5)
#' lam <- ifelse(runif(300) < 0.3, rgamma(300, 2, 1), rgamma(300, 2, 4))
#' safety_gps_fit(rpois(300, lam * E), E)$par
#' @export
safety_gps_fit <- function(n, E, start = c(1/3, 0.2, 0.1, 2, 4)) {
  nll <- function(th) {
    P <- stats::plogis(th[1]); a1 <- exp(th[2]); b1 <- exp(th[3]); a2 <- exp(th[4]); b2 <- exp(th[5])
    l <- log(P * exp(.gps_logmarg(n, E, a1, b1)) + (1 - P) * exp(.gps_logmarg(n, E, a2, b2)) + 1e-300)
    -sum(l)
  }
  th0 <- c(stats::qlogis(start[1]), log(start[-1]))
  best <- NULL
  for (k in 0:3) {
    s <- if (k == 0) th0 else th0 + stats::rnorm(5, 0, 0.3 * k)
    o <- try(stats::optim(s, nll, method = "BFGS", control = list(maxit = 2000)), silent = TRUE)
    if (!inherits(o, "try-error") && (is.null(best) || o$value < best$value)) best <- o
  }
  th <- best$par
  list(par = c(P = stats::plogis(th[1]), a1 = exp(th[2]), b1 = exp(th[3]), a2 = exp(th[4]), b2 = exp(th[5])),
       loglik = -best$value, convergence = best$convergence)
}

# ---------------------------------------------------------------- phase II

#' Economic predictive-probability threshold PP*
#'
#' Continue Phase II if \eqn{PP[\rho V - (1-\rho)c_3] > m c_2}, so
#' \eqn{PP^* = m c_2 / (\rho V - (1-\rho) c_3)}. A non-positive denominator
#' (negative expected value of a Phase II success) returns `Inf`-signed `NA`
#' flag: the value is returned as computed and `viable` is FALSE.
#'
#' @param V Value of an active drug if developed.
#' @param rho Probability Phase III succeeds given Phase II success.
#' @param c3 Cost of a failed Phase III.
#' @param mc2 Cost of finishing Phase II.
#' @return Data frame with `denominator`, `pp_star`, `viable`.
#' @examples
#' p2_pp_star(500, 0.5, 150, 2)
#' @export
p2_pp_star <- function(V, rho, c3, mc2) {
  den <- rho * V - (1 - rho) * c3
  data.frame(V = V, rho = rho, c3 = c3, mc2 = mc2, denominator = den,
             pp_star = mc2 / den, viable = den > 0)
}

#' Bayesian toxicity stopping boundary
#'
#' Stop at look `n` when the number of toxicities `t` satisfies
#' \eqn{\Pr(q > \bar q \mid t, n) > \kappa} with a Beta(a0, b0) prior.
#'
#' @param looks Patient counts at which the rule is checked.
#' @param qbar Maximum acceptable toxicity rate.
#' @param kappa Posterior probability threshold.
#' @param a0,b0 Beta prior.
#' @return Data frame `n`, `stop_at` (smallest toxicity count that stops; `Inf` if none).
#' @examples
#' tox_bayes_boundary()
#' @export
tox_bayes_boundary <- function(looks = seq(5, 40, 5), qbar = 0.2, kappa = 0.75, a0 = 0.5, b0 = 2.5) {
  b <- vapply(looks, function(n) {
    x <- 0:n
    ok <- which(1 - stats::pbeta(qbar, a0 + x, b0 + n - x) > kappa)
    if (length(ok)) min(x[ok]) else Inf
  }, 0)
  data.frame(n = looks, stop_at = b)
}

#' Exact operating characteristics of a toxicity stopping boundary
#'
#' Dynamic programming over the toxicity count: probability of ever stopping and
#' expected number of patients treated (stopping at a look counts that many patients).
#'
#' @param boundary Data frame with `n` and `stop_at` (as from [tox_bayes_boundary()]).
#' @param p True toxicity rates.
#' @return Data frame `p`, `p_stop`, `expected_n`.
#' @examples
#' tox_oc(tox_bayes_boundary(), c(0.1, 0.2, 0.4))
#' @export
tox_oc <- function(boundary, p = c(0.10, 0.15, 0.20, 0.30, 0.40)) {
  nmax <- max(boundary$n)
  res <- t(vapply(p, function(pp) {
    pm <- 1; cur <- 0; ps <- 0; en <- 0  # pm over counts 0..cur (alive); start with n=0
    probs <- 1; prev <- 0
    for (i in seq_len(nrow(boundary))) {
      k <- boundary$n[i] - prev
      for (j in seq_len(k)) probs <- c(probs * (1 - pp), 0) + c(0, probs * pp)
      prev <- boundary$n[i]
      cnt <- seq_along(probs) - 1
      st <- cnt >= boundary$stop_at[i]
      ps_i <- sum(probs[st]); ps <- ps + ps_i; en <- en + ps_i * boundary$n[i]
      probs[st] <- 0
    }
    en <- en + sum(probs) * nmax
    c(ps, en)
  }, c(0, 0)))
  data.frame(p = p, p_stop = res[, 1], expected_n = res[, 2])
}

#' Pocock-type repeated binomial boundary for toxicity monitoring
#'
#' One-sided repeated exact binomial tests of \eqn{H_0: q \le \bar q} at the same
#' looks, with a common nominal level \eqn{a} at every look: stop at look `n` if
#' \eqn{P(X_n \ge t \mid \bar q) \le a}. The nominal `a` is the largest value on a fine
#' grid whose overall stopping probability at `qbar` does not exceed `alpha`.
#'
#' @param looks Look sizes.
#' @param qbar Boundary rate.
#' @param alpha Target overall stopping probability at `qbar`; the default `NULL`
#'   matches the Bayesian rule ([tox_bayes_boundary()] with default arguments).
#' @return List with `boundary` (data frame `n`, `stop_at`), `nominal`, `alpha_achieved`.
#' @examples
#' tox_pocock_boundary(alpha = 0.05)
#' @export
tox_pocock_boundary <- function(looks = seq(5, 40, 5), qbar = 0.2, alpha = NULL) {
  if (is.null(alpha)) alpha <- tox_oc(tox_bayes_boundary(looks, qbar), qbar)$p_stop
  mk <- function(a) {
    b <- vapply(looks, function(n) {
      x <- 0:n; tail <- 1 - stats::pbinom(x - 1, n, qbar)
      ok <- which(tail <= a)
      if (length(ok)) min(x[ok]) else Inf
    }, 0)
    data.frame(n = looks, stop_at = b)
  }
  grid <- exp(seq(log(1e-6), log(0.5), length.out = 600))
  ov <- vapply(grid, function(a) tox_oc(mk(a), qbar)$p_stop, 0)
  ok <- which(ov <= alpha)
  a <- grid[max(ok)]
  list(boundary = mk(a), nominal = a, alpha_achieved = tox_oc(mk(a), qbar)$p_stop)
}

# ---------------------------------------------------------------- BLRM

#' BLRM posterior interval-probability table
#'
#' Two-parameter logistic model \eqn{\mathrm{logit}\,\pi(d) = \log\alpha + \beta\log(d/d^*)}
#' with \eqn{\log\alpha \sim N(m_a, s_a^2)}, \eqn{\log\beta \sim N(m_b, s_b^2)}, posterior
#' computed on a two-dimensional grid (exact up to grid error) or, with `method = "is"`,
#' by prior importance sampling (uses `seed`).
#'
#' @param doses Dose vector.
#' @param n,y Patients and toxicities at each dose (zeros allowed).
#' @param ref Reference dose.
#' @param prior_mean,prior_sd Prior means and sds of (log alpha, log beta).
#' @param cuts Interval cut points (under < cuts[1] <= target < cuts[2] <= over).
#' @param overdose_max Escalation with overdose control: dose admissible if Pr(over) < this.
#' @param method `"grid"` or `"is"`.
#' @param ngrid Grid points per dimension or number of IS draws.
#' @param seed Random seed (importance sampling).
#' @return Data frame: dose, mean toxicity, p_under, p_target, p_over, admissible.
#' @references Neuenschwander B, Branson M, Gsponer T (2008). Statistics in Medicine 27:2420-2439.
#' @examples
#' dose_blrm_table(c(1, 2, 4, 8, 16), n = c(3, 3, 3, 6, 0), y = c(0, 0, 1, 1, 0), ref = 8)
#' @export
dose_blrm_table <- function(doses = c(1, 2, 4, 8, 16), n = c(3, 3, 3, 6, 0), y = c(0, 0, 1, 1, 0),
                            ref = 8, prior_mean = c(stats::qlogis(0.25), 0), prior_sd = c(2, 1),
                            cuts = c(0.16, 0.33), overdose_max = 0.25,
                            method = c("grid", "is"), ngrid = 400, seed = 1) {
  method <- match.arg(method)
  if (method == "grid") {
    la <- seq(prior_mean[1] - 6 * prior_sd[1], prior_mean[1] + 6 * prior_sd[1], length.out = ngrid)
    lb <- seq(prior_mean[2] - 6 * prior_sd[2], prior_mean[2] + 6 * prior_sd[2], length.out = ngrid)
    G <- expand.grid(la = la, lb = lb)
  } else {
    set.seed(seed)
    G <- data.frame(la = stats::rnorm(ngrid, prior_mean[1], prior_sd[1]),
                    lb = stats::rnorm(ngrid, prior_mean[2], prior_sd[2]))
  }
  eta <- function(d) G$la + exp(G$lb) * log(d / ref)
  ll <- 0
  for (i in seq_along(doses)) if (n[i] > 0) ll <- ll + stats::dbinom(y[i], n[i], stats::plogis(eta(doses[i])), log = TRUE)
  lw <- ll + if (method == "grid") stats::dnorm(G$la, prior_mean[1], prior_sd[1], log = TRUE) +
    stats::dnorm(G$lb, prior_mean[2], prior_sd[2], log = TRUE) else 0
  w <- exp(lw - max(lw)); w <- w / sum(w)
  out <- t(vapply(doses, function(d) {
    pi_d <- stats::plogis(eta(d))
    c(sum(w * pi_d), sum(w * (pi_d < cuts[1])), sum(w * (pi_d >= cuts[1] & pi_d < cuts[2])), sum(w * (pi_d >= cuts[2])))
  }, rep(0, 4)))
  data.frame(dose = doses, mean_tox = out[, 1], p_under = out[, 2], p_target = out[, 3],
             p_over = out[, 4], admissible = out[, 4] < overdose_max)
}

# ---------------------------------------------------------------- ECMO

#' ECMO trial data
#'
#' @return Named list of survivors and sample sizes: Michigan (11/11 vs 0/1),
#'   Harvard stage 1 (9/9 vs 6/10), Harvard stage 2 (19/20 on ECMO), UK (deaths 30/93 vs 54/92).
#' @examples
#' ecmo_data()$harvard1
#' @export
ecmo_data <- function() {
  list(michigan = c(yE = 11, nE = 11, yC = 0, nC = 1),
       harvard1 = c(yE = 9, nE = 9, yC = 6, nC = 10),
       harvard12 = c(yE = 28, nE = 29, yC = 6, nC = 10),
       uk = c(yE = 63, nE = 93, yC = 38, nC = 92))
}

#' Pr(ECMO better) with independent Beta priors and Fisher exact test
#'
#' @param d Vector `c(yE, nE, yC, nC)` of survivors and sizes (as from [ecmo_data()]).
#' @param prior Beta(a, b) prior for each survival probability.
#' @return List: `post_prob` = Pr(pE > pC), `fisher_p` one-sided Fisher p value.
#' @examples
#' ecmo_beta_comparison(ecmo_data()$harvard1)
#' @export
ecmo_beta_comparison <- function(d, prior = c(1, 1)) {
  f <- function(x) stats::dbeta(x, prior[1] + d[["yC"]], prior[2] + d[["nC"]] - d[["yC"]]) *
    stats::pbeta(x, prior[1] + d[["yE"]], prior[2] + d[["nE"]] - d[["yE"]], lower.tail = FALSE)
  pp <- stats::integrate(f, 0, 1, rel.tol = 1e-10, subdivisions = 1000L)$value
  tab <- matrix(c(d[["yE"]], d[["nE"]] - d[["yE"]], d[["yC"]], d[["nC"]] - d[["yC"]]), 2, byrow = TRUE)
  list(post_prob = pp, fisher_p = stats::fisher.test(tab, alternative = "greater")$p.value)
}

#' Kass-Greenhouse style prior-sensitivity table for ECMO
#'
#' Logit model \eqn{\eta = \mathrm{logit}\,p_C}, \eqn{\delta = \mathrm{logit}\,p_E - \eta},
#' \eqn{\delta \sim N(0, \sigma_\delta^2)}, control prior on \eqn{\eta}; posterior
#' \eqn{\Pr(\delta > 0 \mid \text{data})} on a 2-D grid.
#'
#' @param d Data vector `c(yE, nE, yC, nC)`.
#' @param sigma_delta Prior sds of the treatment effect.
#' @param control_prior `c(mean, sd)` for eta.
#' @param ngrid Grid size per dimension.
#' @return Data frame `sigma_delta`, `post_prob`.
#' @examples
#' ecmo_prior_sensitivity(ecmo_data()$harvard1)
#' @export
ecmo_prior_sensitivity <- function(d, sigma_delta = c(0.5, 1, 2, 5), control_prior = c(0, 2), ngrid = 1500) {
  eta <- seq(control_prior[1] - 8 * control_prior[2], control_prior[1] + 8 * control_prior[2], length.out = ngrid)
  res <- vapply(sigma_delta, function(s) {
    dl <- seq(-8 * max(s, 1) - 2, 8 * max(s, 1) + 2, length.out = ngrid)
    G <- expand.grid(eta = eta, dl = dl)
    ll <- stats::dbinom(d[["yC"]], d[["nC"]], stats::plogis(G$eta), log = TRUE) +
      stats::dbinom(d[["yE"]], d[["nE"]], stats::plogis(G$eta + G$dl), log = TRUE) +
      stats::dnorm(G$eta, control_prior[1], control_prior[2], log = TRUE) +
      stats::dnorm(G$dl, 0, s, log = TRUE)
    w <- exp(ll - max(ll)); sum(w[G$dl > 0]) / sum(w)
  }, 0)
  data.frame(sigma_delta = sigma_delta, post_prob = res)
}
