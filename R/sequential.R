# Group sequential calculations beyond the shared core (core-sequential.R):
# power under drift, expected information, naive inflation, evidence tolls.

# Gauss-Legendre nodes and weights on [lo, hi]
gs_nodes <- function(lo, hi, n) {
  i <- seq_len(n)
  J <- matrix(0, n, n)
  b <- i[-n] / sqrt(4 * i[-n]^2 - 1)
  J[cbind(i[-n], i[-n] + 1)] <- b
  J[cbind(i[-n] + 1, i[-n])] <- b
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  x <- e$values[o]
  w <- 2 * e$vectors[1, o]^2
  list(x = 0.5 * (hi - lo) * x + 0.5 * (hi + lo), w = 0.5 * (hi - lo) * w)
}

# Gauss-Legendre rule on [-1, 1] (cached)
.gs_rule_cache <- new.env(parent = emptyenv())
gs_rule <- function(n) {
  key <- as.character(n)
  r <- .gs_rule_cache[[key]]
  if (is.null(r)) { r <- gs_nodes(-1, 1, n); .gs_rule_cache[[key]] <- r }
  r
}

# Validate information fractions: finite, positive, strictly increasing, and
# separated enough for the recursion to resolve the transition density.
gs_check_looks <- function(t) {
  if (!is.numeric(t) || length(t) < 1L || anyNA(t) || any(!is.finite(t)))
    stop("`t` must be a non-empty vector of finite information fractions")
  if (any(t <= 0)) stop("information fractions must be positive")
  if (any(diff(t) <= 0)) stop("information fractions must be strictly increasing")
  if (min(diff(c(0, t))) < 1e-9)
    stop("information fractions closer than 1e-9 are not supported")
  invisible(t)
}

# Panel width at each look: 2.5 standard deviations of the shorter of the
# transition into the look (which sets the width of the boundary layer left by
# the previous truncation) and the transition out of it (the kernel that will
# be integrated against the carried density).
gs_panel_width <- function(t) {
  dd <- diff(c(0, t))
  2.5 * sqrt(pmin(dd, c(dd[-1L], Inf)))
}

# Root of f on an adaptively expanded bracket. f must be monotone (or at least
# change sign once near the root); the bracket grows geometrically on both
# sides until the signs differ.
gs_find_root <- function(f, lo, hi, tol = 1e-10, maxit = 30L) {
  flo <- f(lo); fhi <- f(hi); w <- hi - lo; it <- 0L
  while (is.finite(flo) && is.finite(fhi) && sign(flo) == sign(fhi) && flo != 0) {
    it <- it + 1L
    if (it > maxit) stop("could not bracket the root; check the design inputs")
    lo <- lo - w; hi <- hi + w; w <- 2 * w
    flo <- f(lo); fhi <- f(hi)
  }
  if (flo == 0) return(lo)
  if (fhi == 0) return(hi)
  stats::uniroot(f, c(lo, hi), tol = tol)$root
}

# Panelled Gauss-Legendre nodes on [lo, hi]: `npn` nodes per panel, at least
# `nmin` nodes in all, and panels no wider than `wpan`.
gs_panels <- function(lo, hi, wpan, nmin, npn = 12L) {
  P <- max(1, ceiling((hi - lo) / wpan), ceiling(nmin / npn))
  if (P * npn > 6e6) stop("the requested resolution is too fine; looks are too close together")
  r <- gs_rule(npn)
  half <- (hi - lo) / (2 * P)
  mid <- lo + half * (2 * seq_len(P) - 1)
  list(x = rep(mid, each = npn) + half * rep(r$x, P), w = half * rep(r$w, P))
}

# sum_j g_j dnorm((y - x_j - mu) / sd) / sd for ascending x and y, touching
# only the sources within 8.5 sd of each block of targets (banded product).
gs_kernel_apply <- function(y, x, g, mu, sd) {
  out <- numeric(length(y)); nx <- length(x); cut <- 8.5 * sd; chunk <- 256L
  for (s in seq(1L, length(y), by = chunk)) {
    i <- s:min(length(y), s + chunk - 1L)
    j1 <- max(1L, findInterval(y[i[1L]] - mu - cut, x))
    j2 <- min(nx, findInterval(y[i[length(i)]] - mu + cut, x) + 1L)
    if (j2 < j1) next
    j <- j1:j2
    out[i] <- as.vector(stats::dnorm((outer(y[i], x[j], "-") - mu) / sd) %*% g[j]) / sd
  }
  out
}

# Probability (or its log) that the score crosses B at time tk, starting from
# the carried sub-density `st` (NULL state: the first look).
gs_cross <- function(st, tk, B, drift, upper = TRUE, log = FALSE) {
  d <- tk - st$tprev; sd <- sqrt(d); mu <- drift * d
  if (is.null(st$x)) {
    return(stats::pnorm((B - mu) / sd, lower.tail = !upper, log.p = log))
  }
  z <- (B - st$x - mu) / sd
  if (!log) return(sum(st$g * stats::pnorm(z, lower.tail = !upper)))
  lp <- log(st$g) + stats::pnorm(z, lower.tail = !upper, log.p = TRUE)
  m <- max(lp)
  if (!is.finite(m)) return(-Inf)
  m + log(sum(exp(lp - m)))
}

# Carry the sub-density of the score to time tk, restricted to (L, B). The
# integration range is further limited to 9 standard deviations around the
# unconstrained mean (omitted mass below 1e-18); the stopping boundaries
# themselves are never altered. Returns NULL when no mass remains.
gs_advance <- function(st, tk, L, B, drift, wpan, nnode) {
  sk <- sqrt(tk)
  Li <- max(L, drift * tk - 9 * sk); Bi <- min(B, drift * tk + 9 * sk)
  if (!(Li < Bi)) return(NULL)
  nd <- gs_panels(Li, Bi, wpan, nnode)
  d <- tk - st$tprev
  f <- if (is.null(st$x)) stats::dnorm((nd$x - drift * d) / sqrt(d)) / sqrt(d)
       else gs_kernel_apply(nd$x, st$x, st$g, drift * d, sqrt(d))
  list(x = nd$x, g = nd$w * f, tprev = tk)
}

#' Per-look crossing probabilities with optional futility boundary
#'
#' Armitage, McPherson and Rowe (1969) recursion. The score \eqn{S(t)} is
#' Brownian motion with drift `drift` per unit of information fraction; the
#' statistic at look k is \eqn{Z_k = S(t_k)/\sqrt{t_k}}. Unlike
#' [gs_crossing_prob()] it accepts a lower (futility) boundary and reports the
#' mass leaving through it.
#'
#' The sub-density of the score that has not yet stopped is carried between
#' looks on panelled Gauss-Legendre nodes. The panel width at each look is tied
#' to the neighbouring information increments (a panel is at most 2.5 standard
#' deviations of the shorter adjacent transition), so closely spaced or
#' irregular looks are resolved, and `nnode` is a lower bound on the number of nodes per look. The
#' supplied boundaries are applied exactly: the integration range is truncated
#' separately, 9 standard deviations around the unconstrained mean (omitted
#' mass below 1e-18), and a boundary beyond that range changes nothing. A
#' warning is issued if the returned probabilities sum to more than one.
#'
#' @param t Information fractions (strictly increasing, separated by at least 1e-9).
#' @param bounds Upper z boundaries, one per look.
#' @param lower Lower (futility) z boundaries, or `NULL` for none. `-Inf`
#'   entries mean no futility stop at that look. Must not exceed `bounds`.
#' @param drift Expected value of the final z statistic (`theta * sqrt(I_max)`).
#' @param nnode Minimum number of quadrature nodes per look.
#' @return A list with `up` and `down`, per-look probabilities of first
#'   crossing the upper and lower boundary.
#' @examples
#' gs_crossing_two_sided((1:4) / 4, gs_boundary(4), drift = 3)$up
#' @export
gs_crossing_two_sided <- function(t, bounds, lower = NULL, drift = 0, nnode = 160) {
  gs_check_looks(t)
  K <- length(t)
  if (length(bounds) == 1L) bounds <- rep(bounds, K)
  if (length(bounds) != K || anyNA(bounds)) stop("`bounds` must have one non-missing value per look")
  lo <- if (is.null(lower)) rep(-Inf, K) else if (length(lower) == 1L) rep(lower, K) else lower
  if (length(lo) != K || anyNA(lo)) stop("`lower` must be NULL or have one non-missing value per look")
  if (any(lo > bounds)) stop("`lower` boundaries must not exceed the upper boundaries")
  if (length(drift) != 1L || !is.finite(drift)) stop("`drift` must be a single finite number")
  up <- numeric(K); dn <- numeric(K)
  wp <- gs_panel_width(t)
  st <- list(x = NULL, g = NULL, tprev = 0)
  for (k in seq_len(K)) {
    sk <- sqrt(t[k])
    up[k] <- gs_cross(st, t[k], bounds[k] * sk, drift, TRUE)
    if (is.finite(lo[k])) dn[k] <- gs_cross(st, t[k], lo[k] * sk, drift, FALSE)
    if (k < K) {
      st <- gs_advance(st, t[k], lo[k] * sk, bounds[k] * sk, drift, wp[k], nnode)
      if (is.null(st)) break
    }
  }
  if (sum(up) + sum(dn) > 1 + 1e-6)
    warning("crossing probabilities sum to more than one; the quadrature did not converge")
  list(up = pmin(pmax(up, 0), 1), down = pmin(pmax(dn, 0), 1))
}

#' Power of a group sequential design at a given drift
#'
#' @inheritParams gs_crossing_two_sided
#' @return Probability of crossing the upper boundary at some look.
#' @examples
#' gs_power((1:4) / 4, gs_boundary(4), drift = 3.28)
#' @export
gs_power <- function(t, bounds, lower = NULL, drift = 0) {
  sum(gs_crossing_two_sided(t, bounds, lower, drift)$up)
}

#' Drift giving a target power
#'
#' @inheritParams gs_crossing_two_sided
#' @param target Target power (default 0.9), strictly between 0 and 1. The
#'   root is bracketed adaptively, so any attainable target works.
#' @return The drift \eqn{E[Z_K]} at which the overall probability of crossing
#'   the upper boundary equals `target`. The inflation factor over a fixed
#'   design is `(drift / (qnorm(0.975) + qnorm(target)))^2`.
#' @examples
#' gs_power_drift((1:4) / 4, gs_boundary(4))
#' @export
gs_power_drift <- function(t, bounds, lower = NULL, target = 0.9) {
  stopifnot(length(target) == 1L, is.finite(target), target > 0, target < 1)
  gs_find_root(function(d) gs_power(t, bounds, lower, d) - target, 0.5, 8)
}

#' Expected stopping information fraction
#'
#' @inheritParams gs_crossing_two_sided
#' @return \eqn{E[t_\tau]} with early stopping at either boundary.
#' @examples
#' gs_expected_information((1:4) / 4, gs_boundary(4), drift = 3.28)
#' @export
gs_expected_information <- function(t, bounds, lower = NULL, drift = 0) {
  r <- gs_crossing_two_sided(t, bounds, lower, drift)
  stop <- r$up + r$down
  K <- length(t)
  stop[K] <- 1 - sum(stop[-K])
  sum(stop * t)
}

#' Type I error inflation from repeated testing at a fixed cutoff
#'
#' Exact probability that a z statistic exceeds the unadjusted cutoff `z` at
#' least once over `K` equally spaced looks under the null.
#'
#' @param K Number of looks.
#' @param z Unadjusted cutoff (default 1.96).
#' @param t Information fractions (default equally spaced).
#' @return Overall one-sided false positive probability.
#' @examples
#' gs_naive_inflation(2)
#' @export
gs_naive_inflation <- function(K = 4, z = 1.96, t = seq_len(K) / K) {
  if (length(t) == 1L) return(stats::pnorm(z, lower.tail = FALSE))
  gs_crossing_prob(rep(z, length(t)), t)$total
}

#' The three evidence tolls on the |Z| scale
#'
#' Fixed level bar \eqn{z_\alpha}, a Bayes bar with Ockham penalty
#' \eqn{L_n = \log(1 + q)}, \eqn{q = I_n\tau^2}, and an anytime-valid
#' (e-value) bar at threshold \eqn{1/\alpha}. Each of the last two is returned
#' in two clearly labelled forms.
#'
#' The exact finite-information bars come from the normal mixture Bayes factor
#' \eqn{BF_{10} \ge b}, that is
#' \eqn{z^2 \ge (1 + 1/q)(\log(1+q) + 2\log b)}: `bayes_bar_exact` uses
#' \eqn{b = 1} and `anytime_bar_exact` uses \eqn{b = 1/\alpha} (the e-value
#' bar). The asymptotic forms `bayes_bar` \eqn{= \sqrt{L_n}} and `anytime_bar`
#' \eqn{= \sqrt{L_n + 2\log(1/\alpha)}} drop the factor \eqn{1 + 1/q} and
#' agree with the exact bars only when \eqn{q} is large; at small \eqn{q} the
#' asymptotic Bayes bar is far too low (about 0.25 against 1.02 at
#' `info = 1`, `tau = 0.25`).
#'
#' @param info Vector of information levels \eqn{I_n}.
#' @param tau Prior standard deviation.
#' @param alpha One-sided level (default 0.025).
#' @return A data.frame with `info`, `ockham`, `fixed_bar`, the asymptotic
#'   `bayes_bar` and `anytime_bar`, `ratio` (asymptotic anytime over Bayes in
#'   squared units), `surcharge` (\eqn{2\log(1/\alpha)}), and the exact
#'   finite-information `bayes_bar_exact`, `anytime_bar_exact` and
#'   `ratio_exact` (their squared ratio).
#' @references Polson and Sokolov, evidence tolls.
#' @examples
#' gs_tolls(c(25, 100, 1000), tau = 0.25)
#' @export
gs_tolls <- function(info, tau = 0.25, alpha = 0.025) {
  stopifnot(is.numeric(info), all(info > 0), length(tau) == 1L, tau > 0,
            length(alpha) == 1L, alpha > 0, alpha < 1)
  sur <- 2 * log(1 / alpha)
  q <- info * tau^2
  L <- log1p(q)
  ex_b <- bayes_bar(q, 1)
  ex_a <- bayes_bar(q, 1 / alpha)
  data.frame(info = info, ockham = L, fixed_bar = stats::qnorm(1 - alpha),
             bayes_bar = sqrt(L), anytime_bar = sqrt(L + sur),
             ratio = (L + sur) / L, surcharge = sur,
             bayes_bar_exact = sqrt(ex_b), anytime_bar_exact = sqrt(ex_a),
             ratio_exact = ex_a / ex_b)
}
