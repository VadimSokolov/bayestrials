#' Boundary crossing probabilities for a group sequential design
#'
#' Probability that a Brownian motion with drift crosses the upper boundaries
#' `bounds` (z scale) at information fractions `t`, computed with the
#' Armitage, McPherson and Rowe recursion (see [gs_crossing_two_sided()] for
#' the numerical method and its accuracy control).
#'
#' @param bounds Upper z boundaries, one per look.
#' @param t Information fractions (increasing, last equal to 1).
#' @param drift Expected value of the final z statistic (0 under the null).
#' @return A list: `cross` (probability of first crossing at each look) and
#'   `total` (overall crossing probability).
#' @export
gs_crossing_prob <- function(bounds, t, drift = 0) {
  # Armitage, McPherson and Rowe recursion (panelled Gauss-Legendre quadrature
  # whose panel width follows the smallest information increment), exact to
  # quadrature error for irregular and nearly coincident looks.
  r <- gs_crossing_two_sided(t = t, bounds = bounds, lower = NULL, drift = drift, nnode = 200)
  list(cross = r$up, total = sum(r$up))
}

#' Group sequential efficacy boundaries
#'
#' Classical O'Brien-Fleming (`"obf"`, \eqn{c\sqrt{K/k}} for equal spacing,
#' \eqn{c/\sqrt{t_k}} in general) or Pocock (`"pocock"`, constant) boundaries
#' with overall one-sided level `alpha`, solved exactly.
#'
#' @param K Number of looks.
#' @param alpha One-sided type I error (default 0.025).
#' @param type `"obf"` or `"pocock"`.
#' @param t Information fractions (default equally spaced).
#' @return Numeric vector of z boundaries. For `K = 4`, `type = "obf"` this
#'   is 4.049, 2.863, 2.337, 2.024.
#' @export
gs_boundary <- function(K, alpha = 0.025, type = c("obf", "pocock"), t = seq_len(K) / K) {
  type <- match.arg(type)
  gs_check_looks(t)
  stopifnot(length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  shape <- if (type == "obf") 1 / sqrt(t) else rep(1, length(t))
  f <- function(c) gs_crossing_prob(c * shape, t)$total - alpha
  c <- gs_find_root(f, 1, 6)
  c * shape
}

#' Lan-DeMets alpha-spending boundaries
#'
#' Boundaries are found look by look: the sub-density of the Brownian score
#' below the earlier boundaries is carried forward with the same panelled
#' Gauss-Legendre recursion as [gs_crossing_two_sided()], and each boundary is
#' the root of the incremental crossing probability. The spending increments
#' and the root equation are evaluated on the log scale, so looks with
#' negligible spending (for example `t = 0.001` under O'Brien-Fleming spending)
#' give large finite boundaries rather than overflow.
#'
#' @param t Information fractions (strictly increasing, in (0, 1]).
#' @param alpha One-sided level.
#' @param h Retained for backward compatibility and ignored: the recursion no
#'   longer uses a fixed grid step (accuracy is controlled by the panel width).
#' @param spending `"ld_obf"` (O'Brien-Fleming type,
#'   \eqn{2 - 2\Phi(z_{\alpha/2}/\sqrt{t})}) or `"ld_pocock"`
#'   (\eqn{\alpha \log(1 + (e-1)t)}).
#' @return Numeric vector of z boundaries.
#' @export
gs_spending_boundary <- function(t, alpha = 0.025, spending = c("ld_obf", "ld_pocock"),
                                 h = 0.005) {
  spending <- match.arg(spending)
  gs_check_looks(t)
  stopifnot(length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1)
  K <- length(t)
  logcum <- if (spending == "ld_obf") {
    log(2) + stats::pnorm(stats::qnorm(alpha / 2, lower.tail = FALSE) / sqrt(t),
                          lower.tail = FALSE, log.p = TRUE)
  } else {
    log(alpha) + log(log1p((exp(1) - 1) * t))
  }
  gap <- diff(logcum)
  if (any(!is.finite(logcum)) || any(gap <= 0))
    stop("the spending function must be strictly increasing at the supplied looks")
  log1mexp <- function(x) ifelse(x > log(2), log1p(-exp(-x)), log(-expm1(-x)))
  loginc <- c(logcum[1], logcum[-1] + log1mexp(gap))
  b <- numeric(K)
  b[1] <- stats::qnorm(loginc[1], lower.tail = FALSE, log.p = TRUE)
  wp <- gs_panel_width(t)
  st <- gs_advance(list(x = NULL, g = NULL, tprev = 0), t[1], -Inf, b[1] * sqrt(t[1]), 0, wp[1], 200)
  if (K > 1) for (k in 2:K) {
    if (is.null(st)) { b[k:K] <- NA_real_; warning("no probability mass remains after look ", k - 1); break }
    g <- function(bk) pmax(gs_cross(st, t[k], bk * sqrt(t[k]), 0, TRUE, log = TRUE), -1e300) - loginc[k]
    b[k] <- gs_find_root(g, 1, 6)
    if (k < K) st <- gs_advance(st, t[k], -Inf, b[k] * sqrt(t[k]), 0, wp[k], 200)
  }
  b
}
