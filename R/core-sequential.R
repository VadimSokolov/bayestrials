#' Boundary crossing probabilities for a group sequential design
#'
#' Probability that a Brownian motion with drift crosses the upper boundaries
#' `bounds` (z scale) at information fractions `t`, computed with the
#' (computed by the Armitage, McPherson and Rowe recursion).
#'
#' @param bounds Upper z boundaries, one per look.
#' @param t Information fractions (increasing, last equal to 1).
#' @param drift Expected value of the final z statistic (0 under the null).
#' @return A list: `cross` (probability of first crossing at each look) and
#'   `total` (overall crossing probability).
#' @export
gs_crossing_prob <- function(bounds, t, drift = 0) {
  # Armitage, McPherson and Rowe recursion (Gauss-Legendre quadrature), exact to
  # quadrature error and fast for any number of looks; see gs_crossing_two_sided().
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
  shape <- if (type == "obf") 1 / sqrt(t) else rep(1, K)
  f <- function(c) gs_crossing_prob(c * shape, t)$total - alpha
  c <- stats::uniroot(f, c(1, 6), tol = 1e-10)$root
  c * shape
}

#' Lan-DeMets alpha-spending boundaries
#'
#' @param t Information fractions.
#' @param alpha One-sided level.
#' @param h Grid step for the recursion (smaller is more accurate and slower).
#' @param spending `"ld_obf"` (O'Brien-Fleming type,
#'   \eqn{2 - 2\Phi(z_{\alpha/2}/\sqrt{t})}) or `"ld_pocock"`
#'   (\eqn{\alpha \log(1 + (e-1)t)}).
#' @return Numeric vector of z boundaries.
#' @export
gs_spending_boundary <- function(t, alpha = 0.025, spending = c("ld_obf", "ld_pocock"),
                                 h = 0.005) {
  spending <- match.arg(spending)
  cum <- if (spending == "ld_obf") {
    2 * stats::pnorm(stats::qnorm(alpha / 2, lower.tail = FALSE) / sqrt(t), lower.tail = FALSE)
  } else {
    alpha * log(1 + (exp(1) - 1) * t)
  }
  inc <- diff(c(0, cum))
  K <- length(t); b <- numeric(K)
  # Incremental recursion on the Brownian scale W(t) = Z sqrt(t): carry the
  # sub-density of W below the boundary on a uniform grid from look to look.
  b[1] <- stats::qnorm(inc[1], lower.tail = FALSE)
  lo <- -10 * sqrt(max(t))
  grid <- function(top) { g <- seq(lo, top, by = h); if (top - g[length(g)] > 1e-9) g <- c(g, top); g }
  trap <- function(g) { d <- diff(g); c(d, 0) / 2 + c(0, d) / 2 }
  x <- grid(b[1] * sqrt(t[1]))
  f <- stats::dnorm(x, sd = sqrt(t[1]))
  w <- trap(x)
  if (K > 1) for (k in 2:K) {
    sd <- sqrt(t[k] - t[k - 1]); fw <- f * w
    cross <- function(c) sum(fw * stats::pnorm((c - x) / sd, lower.tail = FALSE))
    if (inc[k] < 1e-14) {
      b[k] <- stats::qnorm(inc[k], lower.tail = FALSE)
    } else {
      ck <- stats::uniroot(function(c) cross(c) - inc[k],
                           c(min(x) + 1, 40 * sqrt(t[k])), tol = 1e-10)$root
      b[k] <- ck / sqrt(t[k])
    }
    xn <- grid(b[k] * sqrt(t[k]))
    f <- as.vector(stats::dnorm(outer(xn, x, "-") / sd) %*% fw) / sd
    x <- xn; w <- trap(x)
  }
  b
}
