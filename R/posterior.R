#' Highest posterior density interval of a Beta distribution
#'
#' Grid search over the lower tail probability for the shortest interval with
#' the requested coverage. Valid for unimodal Beta densities (a, b > 1).
#'
#' @param a,b Beta shape parameters (both greater than 1; a warning is
#'   issued otherwise, since the density is not unimodal).
#' @param level Coverage of the interval, strictly between 0 and 1.
#' @param grid Number of grid points for the lower tail probability.
#' @return Named numeric vector with `lower` and `upper`.
#' @examples
#' hpd_beta(4.25, 15.75)
#' @export
hpd_beta <- function(a, b, level = 0.95, grid = 20001) {
  if (!all(is.finite(c(a, b))) || a <= 0 || b <= 0) stop("`a` and `b` must be positive and finite")
  if (a <= 1 || b <= 1)
    warning("hpd_beta() assumes a unimodal density (a, b > 1); the result is not an HPD interval otherwise")
  if (!(level > 0 && level < 1)) stop("`level` must be strictly between 0 and 1")
  lo_q <- seq(0, 1 - level, length.out = grid)
  lo <- stats::qbeta(lo_q, a, b)
  hi <- stats::qbeta(lo_q + level, a, b)
  k <- which.min(hi - lo)
  c(lower = lo[k], upper = hi[k])
}

#' Simulate optional stopping with repeated looks
#'
#' Reproduces the Chapter 5 experiment. (a) Under the null (p1 = p0 = p_null)
#' the one-sided pooled z test is applied at every look, giving the rejection
#' rate at the final look only and at any look. (b) True rates are drawn from
#' independent flat priors, the flat-prior posterior probability Pr(p1 > p0)
#' is computed by Monte Carlo at each look, and the trial stops at the first
#' look where it exceeds `threshold`. All simulation is vectorized over trials.
#'
#' @param n_trials Number of simulated trials.
#' @param n_look Patients per arm added between looks.
#' @param n_looks Number of looks.
#' @param p_null Common response rate in the null scenario.
#' @param threshold Posterior probability threshold for stopping.
#' @param draws Monte Carlo draws per trial and look for the posterior probability.
#' @param seed Optional seed.
#' @return A list with `final_only`, `any_look`, `n_stopped`,
#'   `mean_post_at_stop`, `frac_truly_better`, `frac_truly_better_final_only`
#'   and `mean_post_final_only`.
#' @examples
#' sim_sequential_looks(n_trials = 100, draws = 50, seed = 1)$any_look
#' @export
sim_sequential_looks <- function(n_trials = 4000, n_look = 50, n_looks = 5,
                                 p_null = 0.30, threshold = 0.975,
                                 draws = 1000, seed = NULL) {
  stopifnot(length(n_trials) == 1L, n_trials >= 1, n_trials == round(n_trials),
            length(n_looks) == 1L, n_looks >= 1, n_looks == round(n_looks),
            length(n_look) == 1L, n_look >= 1, n_look == round(n_look),
            p_null >= 0, p_null <= 1, threshold >= 0, threshold <= 1,
            length(draws) == 1L, draws >= 1)
  if (!is.null(seed)) set.seed(seed)
  ns <- n_look * seq_len(n_looks)
  # row-wise cumulative sums that keep the n_trials x n_looks shape even when
  # either dimension is one
  rowcum <- function(m) {
    for (j in seq_len(ncol(m))[-1L]) m[, j] <- m[, j] + m[, j - 1L]
    m
  }
  looks_data <- function(p1, p0) {
    inc1 <- matrix(stats::rbinom(n_trials * n_looks, n_look, p1), n_trials, n_looks)
    inc0 <- matrix(stats::rbinom(n_trials * n_looks, n_look, p0), n_trials, n_looks)
    list(x1 = rowcum(inc1), x0 = rowcum(inc0))
  }
  nmat <- matrix(ns, n_trials, n_looks, byrow = TRUE)
  out <- list()
  # (a) frequentist, null true
  d <- looks_data(rep(p_null, n_trials), rep(p_null, n_trials))
  pbar <- (d$x1 + d$x0) / (2 * nmat)
  se <- sqrt(pmax(pbar * (1 - pbar) * 2 / nmat, 1e-12))
  z <- (d$x1 - d$x0) / nmat / se
  crit <- stats::qnorm(0.975)
  out$final_only <- mean(z[, n_looks] > crit)
  out$any_look <- mean(rowSums(z > crit) > 0)
  # (b) Bayesian, truth from the prior
  p1 <- stats::rbeta(n_trials, 1, 1)
  p0 <- stats::rbeta(n_trials, 1, 1)
  d <- looks_data(p1, p0)
  pr <- matrix(0, n_trials, n_looks)
  for (k in seq_len(n_looks)) {
    a1 <- rep(1 + d$x1[, k], draws); b1 <- rep(1 + ns[k] - d$x1[, k], draws)
    a0 <- rep(1 + d$x0[, k], draws); b0 <- rep(1 + ns[k] - d$x0[, k], draws)
    s1 <- matrix(stats::rbeta(n_trials * draws, a1, b1), n_trials, draws)
    s0 <- matrix(stats::rbeta(n_trials * draws, a0, b0), n_trials, draws)
    pr[, k] <- rowMeans(s1 > s0)
  }
  hit <- pr > threshold
  stopped <- rowSums(hit) > 0
  first <- max.col(hit, ties.method = "first")
  pr_at_stop <- pr[cbind(seq_len(n_trials), first)]
  truth <- p1 > p0
  fin <- pr[, n_looks] > threshold
  out$n_stopped <- sum(stopped)
  out$mean_post_at_stop <- mean(pr_at_stop[stopped])
  out$frac_truly_better <- mean(truth[stopped])
  out$frac_truly_better_final_only <- mean(truth[fin])
  out$mean_post_final_only <- mean(pr[fin, n_looks])
  out
}

#' Beta-binomial probability mass function
#'
#' Predictive pmf of the number of successes in `n` future trials when the
#' success probability has a Beta(`a`, `b`) distribution.
#'
#' @param x Number of successes (vectorized).
#' @param n Number of trials.
#' @param a,b Beta parameters.
#' @param log Return the log pmf?
#' @return The pmf (or log pmf) at `x`; zero outside 0..n.
#' @examples dbetabinom(0:10, 10, 2, 3)
#' @export
dbetabinom <- function(x, n, a, b, log = FALSE) {
  ok <- x >= 0 & x <= n & x == round(x)
  xx <- ifelse(ok, x, 0)
  l <- lchoose(n, xx) + lbeta(xx + a, n - xx + b) - lbeta(a, b)
  l[!ok] <- -Inf
  if (log) l else exp(l)
}

#' Bayes' rule on a grid: prior times likelihood equals posterior
#'
#' Computes the discrete version of Bayes' rule on a grid of parameter values:
#' the posterior is proportional to the prior times the likelihood, normalized
#' to sum to one. The likelihood column is also normalized to sum to one so
#' that the three columns can be drawn on comparable scales.
#'
#' @param theta Grid of parameter values.
#' @param prior Prior weights on the grid (need not sum to one).
#' @param likelihood Likelihood evaluated on the grid, or a function of theta.
#' @return A data.frame with columns `theta`, `prior`, `likelihood` and
#'   `posterior`, each of the last three summing to one.
#' @examples
#' theta <- seq(0, 1, by = 0.1)
#' grid_update(theta, dbeta(theta, 2, 2), dbinom(7, 10, theta))
#' @export
grid_update <- function(theta, prior, likelihood) {
  if (is.function(likelihood)) likelihood <- likelihood(theta)
  stopifnot(length(prior) == length(theta), length(likelihood) == length(theta))
  prior <- prior / sum(prior)
  lik <- likelihood / sum(likelihood)
  post <- prior * likelihood
  if (sum(post) == 0) stop("prior and likelihood have no overlap on the grid")
  data.frame(theta = theta, prior = prior, likelihood = lik, posterior = post / sum(post))
}
