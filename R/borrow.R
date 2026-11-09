# Borrowing historical information: exact operating characteristics for a binary control rate.

.trapz <- function(y, x) sum((y[-1] + y[-length(y)]) / 2 * diff(x))

#' Posterior mixture for a control rate after observing current control data
#'
#' @param yc Number of control responders.
#' @param nc Number of control patients.
#' @param comps Data frame with columns `w`, `a`, `b` (prior mixture of betas).
#' @return Data frame with columns `w`, `a`, `b` for the posterior mixture.
#' @examples
#' borrow_post_control(9, 30, data.frame(w = c(.8, .2), a = c(31, 1), b = c(71, 1)))
#' @export
borrow_post_control <- function(yc, nc, comps) {
  logw <- log(comps$w) + lbeta(comps$a + yc, comps$b + nc - yc) - lbeta(comps$a, comps$b)
  logw <- logw - max(logw)
  pw <- exp(logw)
  pw <- pw / sum(pw)
  data.frame(w = pw, a = comps$a + yc, b = comps$b + nc - yc)
}

#' Posterior probability that treatment beats control under a control mixture
#'
#' Treatment prior is Beta(1, 1). The control posterior is a beta mixture, so
#' \eqn{\Pr(p_t > p_c \mid data) = \sum_i w_i \Pr(p_t > p_{c,i})}, where each term is the
#' exact rank-scale integral of [prob_superior()] (relative tolerance 1e-10). The mixture
#' weights are normalized to sum to one. Vectorized over `yt`.
#'
#' @param yt Treatment responders (vector).
#' @param nt Treatment patients.
#' @param comps_post Data frame with columns `w`, `a`, `b`.
#' @return Numeric vector of \eqn{\Pr(p_t > p_c \mid data)}.
#' @examples
#' borrow_prob_superior(20, 60, data.frame(w = 1, a = 10, b = 22))
#' @export
borrow_prob_superior <- function(yt, nt, comps_post) {
  w <- comps_post$w / sum(comps_post$w)
  out <- numeric(length(yt))
  for (i in seq_along(w)) {
    if (w[i] == 0) next
    out <- out + w[i] * prob_superior(1 + yt, 1 + nt - yt,
                                      comps_post$a[i], comps_post$b[i], method = "exact")
  }
  out
}

#' Control priors compared in the borrowing chapter
#'
#' @param yh,nh Historical control responders and patients.
#' @param a0 Power prior discount.
#' @param w_robust Weight on the informative component of the robust mixture.
#' @return Named list of data frames (`w`, `a`, `b`): No borrowing, Full pooling, Power prior, Robust MAP.
#' @examples
#' borrow_priors(30, 100)
#' @export
borrow_priors <- function(yh, nh, a0 = 0.3, w_robust = 0.8) {
  mk <- function(w, a, b) data.frame(w = w, a = a, b = b)
  out <- list(
    mk(1, 1, 1),
    mk(1, 1 + yh, 1 + nh - yh),
    mk(1, 1 + a0 * yh, 1 + a0 * (nh - yh)),
    mk(c(w_robust, 1 - w_robust), c(1 + yh, 1), c(1 + nh - yh, 1))
  )
  names(out) <- c("No borrowing", "Full pooling", sprintf("Power prior (a0=%g)", a0),
                  sprintf("Robust MAP (w=%g)", w_robust))
  out
}

#' Power prior for a binary control rate
#'
#' @param yh,nh Historical responders and patients.
#' @param a0 Discount, between 0 and 1.
#' @return List with beta parameters `a`, `b`.
#' @examples
#' borrow_power_prior(30, 100, 0.3)
#' @export
borrow_power_prior <- function(yh, nh, a0) {
  list(a = 1 + a0 * yh, b = 1 + a0 * (nh - yh))
}

#' Decision table of a borrowing method
#'
#' Entry `[yc + 1, yt + 1]` is TRUE when \eqn{\Pr(p_t > p_c \mid data) > c}.
#'
#' @param comps Prior mixture for control (data frame `w`, `a`, `b`).
#' @param nt,nc Patients on treatment and control.
#' @param c Decision threshold (named `c` as in the book).
#' @return Logical matrix of dimension (nc + 1) by (nt + 1).
#' @examples
#' borrow_decision_table(data.frame(w = 1, a = 1, b = 1), 10, 5, 0.975)
#' @export
borrow_decision_table <- function(comps, nt, nc, c) {
  tab <- matrix(FALSE, nc + 1, nt + 1)
  for (yc in 0:nc) {
    post <- borrow_post_control(yc, nc, comps)
    tab[yc + 1, ] <- borrow_prob_superior(0:nt, nt, post) > c
  }
  tab
}

#' Exact rejection probability of a decision table
#'
#' @param tab Logical matrix from [borrow_decision_table()].
#' @param pc,pt True control and treatment rates.
#' @param nt,nc Patients on treatment and control.
#' @return Rejection probability (vectorized over `pc` and `pt`).
#' @examples
#' tab <- borrow_decision_table(data.frame(w = 1, a = 1, b = 1), 10, 5, 0.975)
#' borrow_reject_prob(tab, 0.3, 0.3, 10, 5)
#' @export
borrow_reject_prob <- function(tab, pc, pt, nt, nc) {
  n <- max(length(pc), length(pt))
  pc <- rep_len(pc, n); pt <- rep_len(pt, n)
  vapply(seq_len(n), function(i) {
    wc <- dbinom(0:nc, nc, pc[i]); wt <- dbinom(0:nt, nt, pt[i])
    as.numeric(wc %*% tab %*% wt)
  }, numeric(1))
}

#' Posterior weight on the informative component of a robust MAP prior
#'
#' @param yc,nc Current control responders and patients (`yc` may be a vector).
#' @param yh,nh Historical responders and patients.
#' @param w Prior weight on the informative component.
#' @return Posterior weight (numeric vector).
#' @examples
#' borrow_robust_weight(9, 30, 30, 100)
#' @export
borrow_robust_weight <- function(yc, nc, yh, nh, w = 0.8) {
  comps <- borrow_priors(yh, nh, w_robust = w)[[4]]
  vapply(yc, function(k) borrow_post_control(k, nc, comps)$w[1], numeric(1))
}

#' Dynamic effective sample size of a robust prior (mixture-weighted proxy)
#'
#' Convention: the effective sample size (ESS) of a Beta(a, b) component is `a + b`, so the informative component Beta(1 + yh, 1 + nh - yh) has ESS `nh + 2` and the vague Beta(1, 1) component has ESS 2. The function returns the posterior-weight-averaged component ESS, `w_post * (nh + 2) + (1 - w_post) * 2`, where `w_post` is the posterior weight on the informative component after observing `yc` of `nc` current controls.
#'
#' This is a simple proxy, not a general mixture ESS: it is not the moment-matching, predictive-consistency (Morita, Thall and Mueller) or expected-local-information-ratio (Neuenschwander et al.) definition. In particular it equals 2 when `w = 0` (only the vague component remains) and `nh + 2` when `w = 1`. An earlier version multiplied the informative weight by `nh + 2` alone and so returned 0 at `w = 0`.
#'
#' @inheritParams borrow_robust_weight
#' @return Numeric vector.
#' @examples
#' borrow_ess_robust(9, 30, 30, 100)
#' @export
borrow_ess_robust <- function(yc, nc, yh, nh, w = 0.8) {
  wp <- borrow_robust_weight(yc, nc, yh, nh, w)
  wp * (nh + 2) + (1 - wp) * 2
}

#' Hierarchical logit-normal basket model on a grid
#'
#' Half-normal(1) prior on tau, N(0, 2^2) on mu (logit scale), Gauss-type quadrature over basket effects. All likelihoods, quadrature sums and grid weights are computed on the log scale (log-sum-exp), so well-defined posteriors never produce `NaN`, even when some grid points have a vanishing integrated likelihood.
#'
#' @param y,n Responders and patients per basket.
#' @param mu_grid,tau_grid Grids for the logit mean and between-basket sd.
#' @param nq Number of quadrature nodes.
#' @return List with `means` (posterior mean response per basket) and `tau_mean`.
#' @examples
#' borrow_basket_posterior(c(12, 12, 3), c(30, 30, 30), nq = 21)
#' @export
borrow_basket_posterior <- function(y, n, mu_grid = NULL, tau_grid = NULL, nq = 81) {
  K <- length(y)
  if (is.null(mu_grid)) mu_grid <- seq(-3, 1.5, length.out = 46)
  if (is.null(tau_grid)) tau_grid <- seq(0.05, 3.0, length.out = 40)
  z <- seq(-6, 6, length.out = nq)
  lwz <- dnorm(z, log = TRUE); lwz <- lwz - .logsumexp(lwz)
  lp_tau <- log(2) + dnorm(tau_grid, 0, 1, log = TRUE)
  lp_mu <- dnorm(mu_grid, 0, 2, log = TRUE)
  logpost <- matrix(0, length(mu_grid), length(tau_grid))
  ratio <- array(0, c(K, length(mu_grid), length(tau_grid)))
  for (i in seq_along(mu_grid)) {
    for (j in seq_along(tau_grid)) {
      eta <- mu_grid[i] + tau_grid[j] * z
      lth <- plogis(eta, log.p = TRUE); l1th <- plogis(-eta, log.p = TRUE)
      # log binomial likelihood up to the constant choose(n, y), which cancels on normalization
      ll <- outer(y, lth) + outer(n - y, l1th) + rep(lwz, each = K)
      lik <- apply(ll, 1, .logsumexp)
      lnum <- apply(ll + rep(lth, each = K), 1, .logsumexp)
      ratio[, i, j] <- exp(lnum - lik)
      logpost[i, j] <- lp_mu[i] + lp_tau[j] + sum(lik)
    }
  }
  w <- exp(logpost - max(logpost)); w <- w / sum(w)
  means <- vapply(seq_len(K), function(k) sum(w * ratio[k, , ]), numeric(1))
  list(means = means, tau_mean = sum(colSums(w) * tau_grid))
}

.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' Meta-analytic-predictive prior from historical controls (normal-normal, logit scale)
#'
#' Model: empirical logits `est_i = log((y_i + 0.5) / (n_i - y_i + 0.5))` with known variances `se2_i`, `est_i ~ N(theta_i, se2_i)`, `theta_i ~ N(mu, tau^2)`. The prior on the overall logit mean `mu` is flat (improper) and is integrated out analytically, which gives the factor `sqrt(v)` in the marginal posterior of `tau`; the prior on `tau` is half-normal with scale `tau_scale`. The posterior of `tau` is evaluated on a grid with step about 0.005 starting at 0.001.
#'
#' The grid is not fixed: its upper end starts at `max(2, 4 * tau_scale)` and is doubled until the posterior mass beyond it (evaluated by trapezoid on the next doubling) is below `tail_tol` relative to the total, so strongly heterogeneous historical data or a large `tau_scale` are not silently truncated. The retained range and omitted mass are returned as `tau_max` and `tail_mass`.
#'
#' @param y,n Historical responders and patients.
#' @param tau_scale Half-normal scale for the between-trial sd.
#' @param nsamp Number of predictive samples.
#' @param seed Optional seed.
#' @param tail_tol Tolerated omitted posterior mass of `tau`, relative to the total.
#' @return List with `samples` (control rate in a new trial), `beta` (moment-matched `a`, `b`), `taus`, `pw` (posterior grid weights of `tau`), `tau_max` (upper end of the grid) and `tail_mass` (estimated posterior mass of `tau` beyond it).
#' @examples
#' borrow_map_prior(c(24, 33, 27, 45), c(80, 100, 70, 150), nsamp = 1000, seed = 1)$beta
#' @export
borrow_map_prior <- function(y, n, tau_scale = 0.5, nsamp = 200000, seed = NULL, tail_tol = 1e-6) {
  if (!is.null(seed)) set.seed(seed)
  y <- as.numeric(y); n <- as.numeric(n)
  est <- log((y + 0.5) / (n - y + 0.5)); se2 <- 1 / (y + 0.5) + 1 / (n - y + 0.5)
  # log marginal posterior density of tau (flat prior on mu integrated out), mean and variance of mu | tau
  fit <- function(taus) {
    lp <- dnorm(taus, 0, tau_scale, log = TRUE) + log(2)
    mus <- numeric(length(taus)); vmu <- numeric(length(taus))
    for (i in seq_along(taus)) {
      w <- 1 / (se2 + taus[i]^2)
      m <- sum(w * est) / sum(w); v <- 1 / sum(w)
      lp[i] <- lp[i] + 0.5 * log(v) - 0.5 * sum(log(se2 + taus[i]^2)) - 0.5 * sum(w * (est - m)^2)
      mus[i] <- m; vmu[i] <- v
    }
    list(lp = lp, mus = mus, vmu = vmu)
  }
  grid <- function(tmax) seq(0.001, tmax, length.out = max(400L, ceiling(200 * tmax)))
  tmax <- max(2, 4 * tau_scale)
  repeat {
    taus <- grid(tmax)
    f <- fit(taus)
    ext <- seq(tmax, 2 * tmax, length.out = max(400L, ceiling(200 * tmax)))
    lpe <- fit(ext)$lp
    ref <- max(f$lp, lpe)
    tot <- sum(exp(f$lp - ref)) * (taus[2] - taus[1])
    tail_mass <- .trapz(exp(lpe - ref), ext) / (tot + .trapz(exp(lpe - ref), ext))
    if (tail_mass < tail_tol || tmax > 1e4) break
    tmax <- 2 * tmax
  }
  if (tail_mass >= tail_tol) warning("MAP prior: posterior of tau not negligible beyond ", tmax)
  pw <- exp(f$lp - max(f$lp)); pw <- pw / sum(pw)
  idx <- sample.int(length(taus), nsamp, replace = TRUE, prob = pw)
  theta <- rnorm(nsamp, f$mus[idx], sqrt(f$vmu[idx])) + rnorm(nsamp, 0, taus[idx])
  p <- plogis(theta)
  m <- mean(p); v <- var(p)
  k <- m * (1 - m) / v - 1
  list(samples = p, beta = c(a = m * k, b = (1 - m) * k), taus = taus, pw = pw,
       tau_max = tmax, tail_mass = tail_mass)
}
