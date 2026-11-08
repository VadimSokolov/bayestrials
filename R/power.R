#' Squared sum of normal quantiles
#'
#' Internal helper: \eqn{(z_\alpha + z_\beta)^2} with one-sided `alpha`.
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @return Numeric.
#' @noRd
zsum2 <- function(alpha = 0.025, power = 0.90) {
  (stats::qnorm(1 - alpha) + stats::qnorm(power))^2
}

#' Fixed-design sample size from the information equation
#'
#' Solves \eqn{n I \theta^2 = (z_\alpha + z_\beta)^2}.
#'
#' @param theta Effect (for time to event, `-log(HR)`).
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @param info Information per observation `I`.
#' @return Sample size (not rounded).
#' @examples
#' n_fixed(0.2)
#' @export
n_fixed <- function(theta, alpha = 0.025, power = 0.90, info = 1) {
  zsum2(alpha, power) / (info * theta^2)
}

#' Patients per arm for two proportions
#'
#' Unpooled variance, one-sided test.
#'
#' @param p0,p1 Control and treatment response rates.
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @return Patients per arm (not rounded).
#' @examples
#' n_binary(0.3, 0.5)
#' @export
n_binary <- function(p0, p1, alpha = 0.025, power = 0.90) {
  zsum2(alpha, power) * (p0 * (1 - p0) + p1 * (1 - p1)) / (p1 - p0)^2
}

#' Patients per arm for a difference of two means
#'
#' @param delta Difference in means.
#' @param sd Common standard deviation.
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @return Patients per arm (not rounded).
#' @examples
#' n_continuous(5, 10)
#' @export
n_continuous <- function(delta, sd, alpha = 0.025, power = 0.90) {
  2 * zsum2(alpha, power) * sd^2 / delta^2
}

#' Schoenfeld total number of events
#'
#' \eqn{4 (z_\alpha + z_\beta)^2 / \log(HR)^2} for 1:1 allocation.
#'
#' @param hr Hazard ratio.
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @return Total events (not rounded).
#' @references Schoenfeld DA (1983), Biometrics 39:499-503.
#' @examples
#' events_schoenfeld(0.8)
#' @export
events_schoenfeld <- function(hr, alpha = 0.025, power = 0.90) {
  4 * zsum2(alpha, power) / log(hr)^2
}

#' One-sided power of the z test
#'
#' @param theta Effect.
#' @param n Sample size (or events scaled so information is `n * info`).
#' @param alpha One-sided significance level.
#' @param info Information per observation.
#' @return Power.
#' @examples
#' power_z(0.2, 263)
#' @export
power_z <- function(theta, n, alpha = 0.025, info = 1) {
  stats::pnorm(sqrt(n * info) * theta - stats::qnorm(1 - alpha))
}

#' Assurance (Bayesian power)
#'
#' \eqn{P(Z > z_\alpha)} when \eqn{\theta \sim N(m, s^2)} and
#' \eqn{\hat\theta | \theta \sim N(\theta, se^2)}.
#'
#' @param m,s Prior mean and sd of the effect.
#' @param se Standard error of the final estimate.
#' @param alpha One-sided significance level.
#' @return Assurance.
#' @references O'Hagan A, Stevens JW, Campbell MJ (2005), Pharmaceutical Statistics 4:187-201.
#' @examples
#' assurance(0.3, 0.1, 0.15)
#' @export
assurance <- function(m, s, se, alpha = 0.025) {
  stats::pnorm((m - stats::qnorm(1 - alpha) * se) / sqrt(s^2 + se^2))
}

#' Moderate-deviation sample size
#'
#' Solves \eqn{n I \theta^2 = (\sqrt{\log n + C} + z_\beta)^2}, the sample
#' size equation in which the critical value grows with \eqn{\log n}.
#'
#' @param theta Effect.
#' @param C Constant (for example `2 * log(odds)`).
#' @param power Target power.
#' @param info Information per observation.
#' @return Sample size (vectorized over `theta` and `C`).
#' @references Polson N, Sokolov V (2026), working paper on Bayesian standards of evidence.
#' @examples
#' n_md(0.2)
#' @export
n_md <- function(theta, C = 0, power = 0.90, info = 1) {
  zb <- stats::qnorm(power)
  one <- function(theta, C) {
    f <- function(n) n * info * theta^2 - (sqrt(max(log(n) + C, 1e-12)) + zb)^2
    lo <- max(1.0001, exp(-C) + 1e-6)
    stats::uniroot(f, c(lo, 1e12), tol = 1e-12)$root
  }
  d <- data.frame(theta = theta, C = C)
  vapply(seq_len(nrow(d)), function(i) one(d$theta[i], d$C[i]), numeric(1))
}

#' Error rate induced by the moderate-deviation threshold
#'
#' Two-sided error \eqn{2\phi(t)/t} at \eqn{t = \sqrt{\log n + C}}.
#'
#' @param n Sample size.
#' @param C Constant.
#' @return A data.frame with columns `t` and `alpha`.
#' @examples
#' alpha_md(342)
#' @export
alpha_md <- function(n, C = 0) {
  t <- sqrt(log(n) + C)
  data.frame(t = t, alpha = 2 * stats::dnorm(t) / t)
}

#' Moderate-deviation number of events
#'
#' With information `d/4` and \eqn{Z = |\log HR| \sqrt{d}/2}, solves the
#' n-equation on the scale `n = d/4`.
#'
#' @param hr Hazard ratio.
#' @param C Constant.
#' @param power Target power.
#' @return Named numeric vector with `events` and critical value `t`.
#' @examples
#' events_md(0.8)
#' @export
events_md <- function(hr, C = 0, power = 0.90) {
  th <- abs(log(hr)); zb <- stats::qnorm(power)
  f <- function(d) th * sqrt(d) / 2 - (sqrt(max(log(d / 4) + C, 1e-12)) + zb)
  d <- stats::uniroot(f, c(20, 1e8), tol = 1e-12)$root
  c(events = d, t = sqrt(log(d / 4) + C))
}
