#' False positive risk
#'
#' Probability that the drug does not work given a significant result.
#'
#' @param prior Prior probability that the drug works.
#' @param power Power of the test.
#' @param alpha One-sided type I error rate.
#' @return Numeric vector, \eqn{\alpha(1-\pi)/(\alpha(1-\pi) + \beta\pi)}.
#' @references Colquhoun D (2017). The reproducibility of research and the misinterpretation of p-values. Royal Society Open Science.
#' @examples
#' false_positive_risk(0.5, 0.8, 0.025)
#' @export
false_positive_risk <- function(prior, power, alpha = 0.025) {
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
#' @param log_bf Log BF10.
#' @param prior_null Prior probability of the null.
#' @return Posterior probability of the null.
#' @examples
#' posterior_null(1.3, 0.5)
#' @export
posterior_null <- function(log_bf, prior_null = 0.5) {
  odds <- (prior_null / (1 - prior_null)) * exp(-log_bf)
  odds / (1 + odds)
}

#' Lower bound on the Bayes factor for the null
#'
#' `type = "pvalue"`: Sellke-Bayarri-Berger bound `-e p log p` over symmetric
#' unimodal alternatives (valid for p < 1/e). `type = "gaussian"`: bound
#' `exp(-z^2/2)` over all point alternatives, with `x` a z statistic.
#'
#' @param x A p-value (`type = "pvalue"`) or a z statistic (`type = "gaussian"`).
#' @param type Which bound.
#' @return Lower bound on BF01.
#' @references Sellke T, Bayarri MJ, Berger JO (2001). Calibration of p values for testing precise null hypotheses. The American Statistician.
#' @examples
#' min_bf(0.01)
#' @export
min_bf <- function(x, type = c("pvalue", "gaussian")) {
  type <- match.arg(type)
  x <- as.numeric(x)
  if (type == "pvalue") -exp(1) * x * log(x) else exp(-0.5 * x^2)
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
#' @param alpha One-sided level.
#' @param c Bayes factor threshold.
#' @return The value k, or NA if there is no crossing in (1e-3, 1e6).
#' @examples
#' crossing_n()
#' @export
crossing_n <- function(alpha = 0.025, c = 1) {
  z2 <- stats::qnorm(1 - alpha)^2
  f <- function(k) bayes_bar(k, c) - z2
  if (sign(f(1e-3)) == sign(f(1e6))) return(NA_real_)
  stats::uniroot(f, c(1e-3, 1e6), tol = 1e-10)$root
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
