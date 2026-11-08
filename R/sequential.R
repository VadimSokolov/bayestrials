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

#' Per-look crossing probabilities with optional futility boundary
#'
#' Armitage, McPherson and Rowe (1969) recursion with Gauss-Legendre
#' quadrature. The score \eqn{S(t)} is Brownian motion with drift `drift` per
#' unit of information fraction; the statistic at look k is
#' \eqn{Z_k = S(t_k)/\sqrt{t_k}}. Unlike [gs_crossing_prob()] it accepts a lower
#' (futility) boundary and reports the mass leaving through it.
#'
#' @param t Information fractions (increasing).
#' @param bounds Upper z boundaries, one per look.
#' @param lower Lower (futility) z boundaries, or `NULL` for none. `-Inf`
#'   entries mean no futility stop at that look.
#' @param drift Expected value of the final z statistic (`theta * sqrt(I_max)`).
#' @param nnode Number of quadrature nodes.
#' @return A list with `up` and `down`, per-look probabilities of first
#'   crossing the upper and lower boundary.
#' @examples
#' gs_crossing_two_sided((1:4) / 4, gs_boundary(4), drift = 3)$up
#' @export
gs_crossing_two_sided <- function(t, bounds, lower = NULL, drift = 0, nnode = 160) {
  K <- length(t)
  lo <- if (is.null(lower)) rep(-Inf, K) else lower
  up <- numeric(K); dn <- numeric(K)
  x <- NULL; f <- NULL; w <- NULL; tprev <- 0
  cross <- function(tk, B, side) {
    d <- tk - tprev; sd <- sqrt(d); mu <- drift * d
    if (is.null(x)) {
      return(if (side == "up") stats::pnorm((B - mu) / sd, lower.tail = FALSE)
             else stats::pnorm((B - mu) / sd))
    }
    z <- (B - x - mu) / sd
    p <- if (side == "up") stats::pnorm(z, lower.tail = FALSE) else stats::pnorm(z)
    sum(w * f * p)
  }
  for (k in seq_len(K)) {
    sk <- sqrt(t[k])
    B <- bounds[k] * sk
    L <- max(lo[k] * sk, -9 * sk)
    up[k] <- cross(t[k], B, "up")
    if (is.finite(lo[k])) dn[k] <- cross(t[k], L, "down")
    if (k < K) {
      d <- t[k] - tprev; sd <- sqrt(d); mu <- drift * d
      nd <- gs_nodes(L, B, nnode)
      if (is.null(x)) {
        fy <- stats::dnorm((nd$x - mu) / sd) / sd
      } else {
        Km <- stats::dnorm((outer(nd$x, x, "-") - mu) / sd) / sd
        fy <- as.vector(Km %*% (w * f))
      }
      x <- nd$x; w <- nd$w; f <- fy; tprev <- t[k]
    }
  }
  list(up = up, down = dn)
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
#' @param target Target power (default 0.9).
#' @return The drift \eqn{E[Z_K]} at which the overall probability of crossing
#'   the upper boundary equals `target`. The inflation factor over a fixed
#'   design is `(drift / (qnorm(0.975) + qnorm(target)))^2`.
#' @examples
#' gs_power_drift((1:4) / 4, gs_boundary(4))
#' @export
gs_power_drift <- function(t, bounds, lower = NULL, target = 0.9) {
  stats::uniroot(function(d) gs_power(t, bounds, lower, d) - target,
                 c(0.5, 8), tol = 1e-10)$root
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
#' Fixed level bar \eqn{z_\alpha}, Bayes bar \eqn{\sqrt{L_n}} with Ockham
#' penalty \eqn{L_n = \log(1 + I_n\tau^2)}, and anytime-valid bar
#' \eqn{\sqrt{L_n + 2\log(1/\alpha)}}.
#'
#' @param info Vector of information levels \eqn{I_n}.
#' @param tau Prior standard deviation.
#' @param alpha One-sided level (default 0.025).
#' @return A data.frame with `info`, `ockham`, `fixed_bar`, `bayes_bar`,
#'   `anytime_bar`, `ratio` (anytime over Bayes in squared units) and
#'   `surcharge`.
#' @references Polson and Sokolov, evidence tolls.
#' @examples
#' gs_tolls(c(25, 100, 1000), tau = 0.25)
#' @export
gs_tolls <- function(info, tau = 0.25, alpha = 0.025) {
  sur <- 2 * log(1 / alpha)
  L <- log(1 + info * tau^2)
  data.frame(info = info, ockham = L, fixed_bar = stats::qnorm(1 - alpha),
             bayes_bar = sqrt(L), anytime_bar = sqrt(L + sur),
             ratio = (L + sur) / L, surcharge = sur)
}
