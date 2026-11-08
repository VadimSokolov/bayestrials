#' Beta posterior for a binomial outcome
#'
#' @param x Number of responders.
#' @param n Number of patients.
#' @param a,b Beta prior parameters (default uniform).
#' @return A list with posterior parameters `a` and `b`.
#' @export
beta_post <- function(x, n, a = 1, b = 1) {
  list(a = a + x, b = b + n - x)
}

#' Posterior probability that one Beta-distributed rate exceeds another
#'
#' Computes \eqn{\Pr(p_1 > p_0 + m)} for independent \eqn{p_1 \sim Beta(a_1, b_1)}
#' and \eqn{p_0 \sim Beta(a_0, b_0)}. Vectorized over all arguments.
#'
#' @param a1,b1 Beta parameters of the treatment arm.
#' @param a0,b0 Beta parameters of the control arm.
#' @param margin Margin `m` (default 0).
#' @param method `"exact"` (numerical integration), `"normal"` (normal
#'   approximation to the difference, fast, used in large simulations) or
#'   `"mc"` (Monte Carlo).
#' @param ndraw Number of Monte Carlo draws for `method = "mc"`.
#' @return Numeric vector of probabilities.
#' @export
prob_superior <- function(a1, b1, a0, b0, margin = 0,
                          method = c("exact", "normal", "mc"), ndraw = 1e5) {
  method <- match.arg(method)
  args <- data.frame(a1 = a1, b1 = b1, a0 = a0, b0 = b0, m = margin)
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
    with(args[i, ], {
      f <- function(p0) stats::dbeta(p0, a0, b0) *
        stats::pbeta(pmin(pmax(p0 + m, 0), 1), a1, b1, lower.tail = FALSE)
      stats::integrate(f, 0, 1, rel.tol = 1e-10, subdivisions = 1000L)$value
    })
  }, numeric(1))
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
