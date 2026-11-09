#' Beta posterior for a binomial outcome
#'
#' @param x Number of responders.
#' @param n Number of patients.
#' @param a,b Beta prior parameters (default uniform), positive and finite.
#' @return A list with posterior parameters `a` and `b`. An error is raised
#'   unless `n` is a nonnegative whole number and `x` a whole number in `[0, n]`.
#' @export
beta_post <- function(x, n, a = 1, b = 1) {
  if (!is.numeric(x) || !is.numeric(n) || anyNA(x) || anyNA(n) ||
      any(x < 0) || any(n < 0) || any(x > n) || any(x != round(x)) || any(n != round(n)))
    stop("`x` and `n` must be whole numbers with 0 <= x <= n")
  if (!is.numeric(a) || !is.numeric(b) || anyNA(a) || anyNA(b) ||
      any(!is.finite(a)) || any(!is.finite(b)) || any(a <= 0) || any(b <= 0))
    stop("prior shapes `a` and `b` must be positive and finite")
  list(a = a + x, b = b + n - x)
}

#' Posterior probability that one Beta-distributed rate exceeds another
#'
#' Computes \eqn{\Pr(p_1 > p_0 + m)} for independent \eqn{p_1 \sim Beta(a_1, b_1)}
#' and \eqn{p_0 \sim Beta(a_0, b_0)}. Vectorized over all arguments.
#'
#' The `"exact"` method integrates over the quantile (rank) scale of the
#' *narrower* of the two distributions, \eqn{\int_0^1 G(Q(u))\,du}, where `Q`
#' is its Beta quantile function and `G` the other arm's Beta (survival)
#' function at the shifted point. On this scale the integrand is bounded,
#' monotone and smooth for any shapes, including U-shaped, concentrated and
#' endpoint-heavy Betas, whereas a fixed grid or an adaptive integral in the
#' probability scale can miss a concentrated peak or endpoint mass. The
#' integration is adaptive with relative tolerance 1e-10; where the
#' integrand has a kink (for example a margin that clips at 0 or 1) it is split
#' there.
#'
#' @param a1,b1 Beta parameters of the treatment arm (positive, finite).
#' @param a0,b0 Beta parameters of the control arm.
#' @param margin Margin `m` (default 0).
#' @param method `"exact"` (numerical integration as above), `"normal"` (normal
#'   approximation to the difference, fast, used in large simulations) or
#'   `"mc"` (Monte Carlo).
#' @param ndraw Number of Monte Carlo draws for `method = "mc"`.
#' @return Numeric vector of probabilities.
#' @examples
#' prob_superior(5, 7, 5, 7)            # exchangeable: exactly 0.5
#' prob_superior(0.01, 0.01, 0.01, 0.01) # U-shaped: still 0.5
#' @export
prob_superior <- function(a1, b1, a0, b0, margin = 0,
                          method = c("exact", "normal", "mc"), ndraw = 1e5) {
  method <- match.arg(method)
  args <- data.frame(a1 = a1, b1 = b1, a0 = a0, b0 = b0, m = margin)
  shp <- as.matrix(args[, c("a1", "b1", "a0", "b0")])
  if (anyNA(args) || any(!is.finite(shp)) || any(shp <= 0))
    stop("Beta shape parameters must be positive and finite")
  if (any(!is.finite(args$m))) stop("`margin` must be finite")
  if (method == "normal") {
    m1 <- args$a1 / (args$a1 + args$b1); m0 <- args$a0 / (args$a0 + args$b0)
    v1 <- m1 * (1 - m1) / (args$a1 + args$b1 + 1)
    v0 <- m0 * (1 - m0) / (args$a0 + args$b0 + 1)
    return(stats::pnorm((m1 - m0 - args$m) / sqrt(v1 + v0)))
  }
  if (method == "mc") {
    return(vapply(seq_len(nrow(args)), function(i) {
      with(args[i, ], mean(stats::rbeta(ndraw, a1, b1) > stats::rbeta(ndraw, a0, b0) + m))
    }, numeric(1)))
  }
  vapply(seq_len(nrow(args)), function(i) {
    with(args[i, ], prob_superior_exact(a1, b1, a0, b0, m))
  }, numeric(1))
}

# Exact Pr(p1 > p0 + m) by integrating over the rank scale of the narrower arm.
# Each half of the rank scale is handled from its own end of [0, 1]: the lower
# half uses the quantile of p, the upper half the quantile of 1 - p (a Beta
# with swapped shapes), so mass within 1e-300 of either endpoint is resolved
# instead of being rounded onto 0 or 1.
prob_superior_exact <- function(a1, b1, a0, b0, m) {
  v <- function(a, b) a * b / ((a + b)^2 * (a + b + 1))
  if (m >= 1) return(0)
  if (m <= -1) return(1)
  # point (p, x = 1 - p) on the narrower arm's quantile scale, each coordinate accurate at its own end
  pt <- function(r, a, b, lower) {
    if (lower) { p <- stats::qbeta(r, a, b); x <- stats::qbeta(r, b, a, lower.tail = FALSE) }
    else { p <- stats::qbeta(r, a, b, lower.tail = FALSE); x <- stats::qbeta(r, b, a) }
    small <- p < 0.5
    list(p = ifelse(small, p, 1 - x), x = ifelse(small, 1 - p, x), small = small)
  }
  # Pr(other > t) (upper) or Pr(other < t) (lower) for t = p + shift, using the accurate end
  tail_at <- function(q, shift, a, b, upper) {
    tp <- ifelse(q$small, q$p + shift, 1 - (q$x - shift))
    tx <- ifelse(q$small, 1 - (q$p + shift), q$x - shift)
    tp <- pmin(pmax(tp, 0), 1); tx <- pmin(pmax(tx, 0), 1)
    useP <- tp < 0.5
    if (upper) ifelse(useP, stats::pbeta(tp, a, b, lower.tail = FALSE), stats::pbeta(tx, b, a))
    else ifelse(useP, stats::pbeta(tp, a, b), stats::pbeta(tx, b, a, lower.tail = FALSE))
  }
  if (v(a0, b0) <= v(a1, b1)) {
    # condition on p0 (narrower): Pr(p1 > p0 + m | p0)
    gl <- function(u) tail_at(pt(u, a0, b0, TRUE), m, a1, b1, TRUE)
    gu <- function(w) tail_at(pt(w, a0, b0, FALSE), m, a1, b1, TRUE)
    kinks <- c(stats::pbeta(max(1 - m, 0), a0, b0), stats::pbeta(max(-m, 0), a0, b0))
  } else {
    # condition on p1 (narrower): Pr(p0 < p1 - m | p1)
    gl <- function(u) tail_at(pt(u, a1, b1, TRUE), -m, a0, b0, FALSE)
    gu <- function(w) tail_at(pt(w, a1, b1, FALSE), -m, a0, b0, FALSE)
    kinks <- c(stats::pbeta(max(m, 0), a1, b1), stats::pbeta(min(1 + m, 1), a1, b1))
  }
  half <- function(g, cuts) {
    cuts <- sort(unique(c(0, 0.5, cuts[cuts > 0 & cuts < 0.5])))
    tot <- 0
    for (k in seq_len(length(cuts) - 1L)) {
      tot <- tot + stats::integrate(g, cuts[k], cuts[k + 1L], rel.tol = 1e-10,
                                    abs.tol = 1e-13, subdivisions = 500L,
                                    stop.on.error = FALSE)$value
    }
    tot
  }
  min(max(half(gl, kinks) + half(gu, 1 - kinks), 0), 1)
}

#' Conjugate normal update
#'
#' Combines a normal prior \eqn{N(\mu_0, \sigma_0^2)} with an estimate
#' \eqn{\hat\theta} of standard error `se` (normal likelihood).
#'
#' @param prior_mean,prior_sd Prior mean and standard deviation.
#' @param est,se Estimate and its standard error.
#' @return A list with posterior `mean` and `sd`.
#' @export
normal_update <- function(prior_mean, prior_sd, est, se) {
  w0 <- 1 / prior_sd^2; w1 <- 1 / se^2
  list(mean = (w0 * prior_mean + w1 * est) / (w0 + w1), sd = sqrt(1 / (w0 + w1)))
}

#' Effective sample size of a prior
#'
#' `ess_beta()` returns \eqn{a + b} (Morita, Thall and Mueller 2008).
#' `ess_normal()` returns \eqn{\sigma^2 / \sigma_0^2}, the number of
#' observations of per-observation standard deviation `sigma` the prior is worth;
#' for a log hazard ratio with 1:1 allocation use `sigma = 2` (events).
#'
#' @param a,b Beta parameters.
#' @param prior_sd Prior standard deviation.
#' @param sigma Per-observation standard deviation.
#' @export
ess_beta <- function(a, b) a + b

#' @rdname ess_beta
#' @export
ess_normal <- function(prior_sd, sigma) sigma^2 / prior_sd^2
