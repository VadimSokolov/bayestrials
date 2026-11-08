#' Effective number of events in a normal prior
#'
#' A prior \eqn{N(., s_0^2)} on the log hazard ratio is worth `4 / s0^2`
#' events under 1:1 allocation.
#'
#' @param s0 Prior standard deviation.
#' @return Events.
#' @examples
#' ess_events(0.5)
#' @export
ess_events <- function(s0) 4 / s0^2

#' Posterior for a mixture of normal priors
#'
#' @param est,se Estimate and standard error.
#' @param comps A data.frame (or list) with columns `weight`, `mean`, `sd`.
#' @return A list with posterior `weights`, `means`, `sds` and `p_pos`,
#'   the posterior probability that theta > 0.
#' @examples
#' mixture_posterior(0.2, 0.1, data.frame(weight = c(.5, .5), mean = c(0, .3), sd = c(.5, .2)))
#' @export
mixture_posterior <- function(est, se, comps) {
  comps <- as.data.frame(comps)
  w <- comps$weight * stats::dnorm(est, comps$mean, sqrt(comps$sd^2 + se^2))
  w <- w / sum(w)
  up <- normal_update(comps$mean, comps$sd, est, se)
  list(weights = w, means = up$mean, sds = up$sd,
       p_pos = sum(w * stats::pnorm(0, up$mean, up$sd, lower.tail = FALSE)))
}

#' Conditional means prior draws for a logistic dose-response
#'
#' Prior for logit p = b0 + slope * (x - 30), set by Beta(a1, b1) at dose
#' `x1` and Beta(a2, b2) at dose `x2`.
#'
#' @param a1,b1 Beta parameters at `x1`.
#' @param a2,b2 Beta parameters at `x2`.
#' @param x1,x2 Doses.
#' @param n Number of draws.
#' @param x0 Centering dose (default 30).
#' @param seed Optional seed.
#' @return A data.frame with columns `b0` and `slope`.
#' @references Bedrick EJ, Christensen R, Johnson W (1996), JASA 91:1450-1460.
#' @examples
#' head(cmp_samples(1, 4, 4, 1, 10, 50, 5, seed = 1), 2)
#' @export
cmp_samples <- function(a1, b1, a2, b2, x1, x2, n, x0 = 30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  l1 <- stats::qlogis(stats::rbeta(n, a1, b1))
  l2 <- stats::qlogis(stats::rbeta(n, a2, b2))
  slope <- (l2 - l1) / (x2 - x1)
  data.frame(b0 = l1 + slope * (x0 - x1), slope = slope)
}

#' Probability of success (assurance) of a z test
#'
#' @param m,s Prior mean and sd of theta (`s = 0` gives ordinary power).
#' @param se Standard error of the final estimate.
#' @param z_crit Critical value.
#' @return Probability that the test rejects.
#' @examples
#' prob_success(0.3, 0.1, 0.2)
#' @export
prob_success <- function(m, s, se, z_crit = 1.959964) {
  stats::pnorm((m - z_crit * se) / sqrt(s^2 + se^2))
}
