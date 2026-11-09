#' False positive risk
#'
#' Probability that the drug does not work given a significant result.
#'
#' @param prior Prior probability that the drug works.
#' @param power Power of the test.
#' @param alpha One-sided type I error rate.
#' @return Numeric vector, \eqn{\alpha(1-\pi)/(\alpha(1-\pi) + (1-\beta)\pi)},
#'   where \eqn{\pi} is `prior`, \eqn{\alpha} the type I error and
#'   \eqn{1-\beta} = `power` (the code takes `power`, not the type II error
#'   rate \eqn{\beta}).
#' @references Colquhoun D (2017). The reproducibility of research and the misinterpretation of p-values. Royal Society Open Science.
#' @examples
#' false_positive_risk(0.5, 0.8, 0.025)
#' @export
false_positive_risk <- function(prior, power, alpha = 0.025) {
  if (anyNA(c(prior, power, alpha)) || any(c(prior, power, alpha) < 0) ||
      any(c(prior, power, alpha) > 1))
    stop("`prior`, `power` and `alpha` must be probabilities in [0, 1]")
  alpha * (1 - prior) / (alpha * (1 - prior) + power * prior)
}

#' Power of the one-sided z test
#'
#' @param theta_over_se True effect in standard errors.
#' @param alpha One-sided level.
#' @return Power.
#' @examples
#' power_one_sided(2.8)
#' @export
power_one_sided <- function(theta_over_se, alpha = 0.025) {
  1 - stats::pnorm(stats::qnorm(1 - alpha) - theta_over_se)
}

#' Bayes factor for a normal prior on the effect against a point null
#'
#' H1: theta ~ N(0, tau^2) against H0: theta = 0, with `k = I_n * tau^2` the
#' prior-to-data information ratio.
#'
#' @param z Z statistic.
#' @param k Information ratio.
#' @param log If TRUE return the log Bayes factor.
#' @return BF10 (or its log).
#' @references Lindley DV (1957). A statistical paradox. Biometrika.
#' @examples
#' bf10_normal(2.5, 10)
#' @export
bf10_normal <- function(z, k, log = FALSE) {
  lb <- -0.5 * log1p(k) + 0.5 * z^2 * k / (1 + k)
  if (log) lb else exp(lb)
}

#' Posterior probability of the null
#'
#' @param log_bf Log BF10 (any real value, evaluated on the log-odds scale so
#'   that very large or small Bayes factors do not overflow).
#' @param prior_null Prior probability of the null.
#' @return Posterior probability of the null.
#' @examples
#' posterior_null(1.3, 0.5)
#' @export
posterior_null <- function(log_bf, prior_null = 0.5) {
  if (anyNA(prior_null) || any(prior_null < 0) || any(prior_null > 1))
    stop("`prior_null` must be a probability in [0, 1]")
  # logistic form of the posterior odds: stable for |log_bf| in the thousands
  stats::plogis(stats::qlogis(prior_null) - log_bf)
}

#' Lower bound on the Bayes factor for the null
#'
#' `type = "pvalue"`: Sellke-Bayarri-Berger bound \eqn{-e\,p\log p} over
#' symmetric unimodal alternatives for \eqn{p < 1/e}, extended by the value 1
#' for \eqn{p \ge 1/e} (where no evidence against the null is possible under
#' the calibration); the limit 0 is returned at \eqn{p = 0}.
#' `type = "gaussian"`: bound `exp(-z^2/2)` over all point alternatives, with
#' `x` a z statistic.
#'
#' @param x A p-value in [0, 1] (`type = "pvalue"`) or a z statistic
#'   (`type = "gaussian"`).
#' @param type Which bound.
#' @return Lower bound on BF01.
#' @references Sellke T, Bayarri MJ, Berger JO (2001). Calibration of p values for testing precise null hypotheses. The American Statistician.
#' @examples
#' min_bf(0.01)
#' @export
min_bf <- function(x, type = c("pvalue", "gaussian")) {
  type <- match.arg(type)
  x <- as.numeric(x)
  if (type == "gaussian") return(exp(-0.5 * x^2))
  if (any(x < 0 | x > 1, na.rm = TRUE)) stop("p-values must lie in [0, 1]")
  out <- rep(1, length(x))
  lt <- !is.na(x) & x < exp(-1)
  out[lt] <- ifelse(x[lt] == 0, 0, -exp(1) * x[lt] * log(x[lt]))
  out[is.na(x)] <- NA_real_
  out
}

#' Squared-z bar of the mixture Bayes rule
#'
#' @param k Information ratio `I_n * tau^2`.
#' @param c Bayes factor threshold.
#' @return \eqn{(1+k)/k (\log(1+k) + 2 \log c)}.
#' @examples
#' bayes_bar(10, 4)
#' @export
bayes_bar <- function(k, c = 1) (1 + k) / k * (log1p(k) + 2 * log(c))

#' Information ratio at which the Bayes bar equals the fixed-level bar
#'
#' Solves \eqn{(1+k)/k\,(\log(1+k) + 2\log c) = z_\alpha^2} for the
#' information ratio `k`. For `c = 1` (and `c < 1`) the Bayes bar increases
#' through the fixed bar once (the Lindley regime), so there is one root. For
#' `c > 1` the bar is large at small `k`, dips, and rises again: there can be
#' two roots, or none if the dip stays above the fixed bar. The scalar value is
#' the **upper** crossover, the `k` beyond which the Bayes bar stays above the
#' fixed bar; `all = TRUE` returns every root.
#'
#' @param alpha One-sided level, in (0, 1).
#' @param c Bayes factor threshold (positive).
#' @param all If `TRUE`, return all roots in increasing order (possibly
#'   `numeric(0)`) instead of the upper one.
#' @return The upper crossover `k`, or `NA` if there is no crossing in
#'   (1e-8, 1e12) (`all = FALSE`); a numeric vector of roots (`all = TRUE`).
#' @examples
#' crossing_n()
#' crossing_n(0.025, 2, all = TRUE)   # two crossings: 1.489 and 5.378
#' @export
crossing_n <- function(alpha = 0.025, c = 1, all = FALSE) {
  stopifnot(length(alpha) == 1L, is.finite(alpha), alpha > 0, alpha < 1,
            length(c) == 1L, is.finite(c), c > 0)
  z2 <- stats::qnorm(1 - alpha)^2
  f <- function(k) bayes_bar(k, c) - z2
  lk <- seq(log(1e-8), log(1e12), length.out = 2001)
  # the bar can dip between grid points: add the minimizer of f on the log scale
  fl <- function(u) f(exp(u))
  imin <- which.min(fl(lk))
  if (imin > 1 && imin < length(lk)) {
    o <- stats::optimize(fl, c(lk[imin - 1], lk[imin + 1]), tol = 1e-12)
    lk <- sort(c(lk, o$minimum))
  }
  fv <- fl(lk)
  roots <- numeric(0)
  for (i in which(fv[-length(fv)] * fv[-1] < 0 | fv[-length(fv)] == 0)) {
    roots <- c(roots, if (fv[i] == 0) exp(lk[i]) else
      exp(stats::uniroot(fl, lk[c(i, i + 1)], tol = 1e-13)$root))
  }
  roots <- sort(unique(roots))
  if (all) return(roots)
  if (length(roots)) max(roots) else NA_real_
}

#' Estimate threshold for the flat-prior success criterion
#'
#' With a flat prior, Pr(theta > a | data) > c iff the estimate exceeds
#' `a + qnorm(c) * se`.
#'
#' @param a Threshold effect.
#' @param c Posterior probability cutoff.
#' @param se Standard error.
#' @return Estimate threshold.
#' @examples
#' threshold_flat(0.1, 0.95, 0.05)
#' @export
threshold_flat <- function(a, c, se) a + stats::qnorm(c) * se

#' Type I error of a posterior-probability criterion
#'
#' One-sided type I error at theta = 0 for Pr(theta > a | data) > c, under a
#' flat prior (`type1_posterior_flat`) or a N(0, tau^2) prior
#' (`type1_posterior_normal`).
#'
#' @param a Threshold effect.
#' @param c Posterior probability cutoff.
#' @param se Standard error of the estimate.
#' @param tau Prior standard deviation (normal prior).
#' @return Type I error.
#' @examples
#' type1_posterior_flat(0.1, 0.95, 0.05)
#' @export
type1_posterior_flat <- function(a, c, se) {
  1 - stats::pnorm(threshold_flat(a, c, se) / se)
}

#' @rdname type1_posterior_flat
#' @export
type1_posterior_normal <- function(a, c, se, tau) {
  w <- tau^2 / (tau^2 + se^2)
  v <- w * se^2
  thr <- (a + stats::qnorm(c) * sqrt(v)) / w
  1 - stats::pnorm(thr / se)
}

#' Power of the posterior-probability criterion under a normal prior
#'
#' @param a Threshold effect.
#' @param c Posterior probability cutoff.
#' @param se Standard error.
#' @param tau Prior standard deviation.
#' @param theta True effect.
#' @return Power.
#' @examples
#' power_posterior_normal(0.1, 0.95, 0.05, 0.3, 0.15)
#' @export
power_posterior_normal <- function(a, c, se, tau, theta) {
  w <- tau^2 / (tau^2 + se^2)
  v <- w * se^2
  thr <- (a + stats::qnorm(c) * sqrt(v)) / w
  1 - stats::pnorm((thr - theta) / se)
}

#' Cutoff c giving type I error alpha at theta = 0
#'
#' `calibrate_c_flat()` for a flat prior, `calibrate_c_normal()` for a
#' N(0, tau^2) prior.
#'
#' @param a Threshold effect.
#' @param se Standard error.
#' @param tau Prior standard deviation.
#' @param alpha One-sided type I error.
#' @return The posterior probability cutoff c.
#' @examples
#' calibrate_c_flat(0.1, 0.05)
#' @export
calibrate_c_flat <- function(a, se, alpha = 0.025) {
  stats::pnorm(stats::qnorm(1 - alpha) - a / se)
}

#' @rdname calibrate_c_flat
#' @export
calibrate_c_normal <- function(a, se, tau, alpha = 0.025) {
  w <- tau^2 / (tau^2 + se^2)
  v <- w * se^2
  thr <- stats::qnorm(1 - alpha) * se
  stats::pnorm((w * thr - a) / sqrt(v))
}
