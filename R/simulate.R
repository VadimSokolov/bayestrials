# Operating-characteristics simulator for a two-arm binary trial with interim looks.
# Port of code/ch12_sim.py. All simulation is matrix-vectorized over trials: a simulated trial
# costs two binomial draws per look and one table lookup.

.oc_cache <- new.env(parent = emptyenv())

#' Posterior probability table for a two-arm binary trial
#'
#' Returns the matrix `T[x0 + 1, x1 + 1] = Pr(p1 - p0 > a | x0, x1)` for `x0`, `x1` successes
#' out of `n` per arm with independent Beta(`prior`) priors on each arm. Each entry is the exact
#' rank-scale integral of [prob_superior()] (relative tolerance 1e-10), so exchangeable states
#' give exactly 0.5 and concentrated or U-shaped posteriors are not truncated. When `a = 0`
#' the table is exchangeable and only half of it is integrated. It is computed once per look, so
#' that a simulated trial reads its posterior probability by table lookup. Tables are cached
#' within the session under a key that records `n`, `a` and `prior` exactly (to full double
#' precision), so distinct arguments never share a table.
#'
#' Cost grows like `n^2` integrals (about 0.6 ms each; a table for `n = 160` takes about 5 s
#' with `a = 0` and 15 s otherwise). On Unix, `options(bayestrials.cores = 4)` computes the rows
#' of a table in parallel with forked workers.
#'
#' @param n Patients per arm at the look.
#' @param a Superiority margin (default 0).
#' @param prior Beta prior `c(shape1, shape2)` used on both arms.
#' @param G Deprecated and ignored (earlier versions used a midpoint grid of `G` points).
#' @return Numeric matrix of dimension `(n + 1) x (n + 1)`, rows indexed by control
#'   successes and columns by treatment successes (offset by one).
#' @export
#' @examples
#' oc_post_table(10)[6, 9]
oc_post_table <- function(n, a = 0, prior = c(1, 1), G = NULL) {
  if (length(n) != 1L || !is.numeric(n) || !is.finite(n) || n < 0 || n != round(n))
    stop("`n` must be a single nonnegative whole number")
  if (length(prior) != 2L || !is.numeric(prior) || any(!is.finite(prior)) || any(prior <= 0))
    stop("`prior` must be two positive, finite Beta shapes")
  if (length(a) != 1L || !is.numeric(a) || !is.finite(a)) stop("`a` must be a single finite number")
  key <- paste(sprintf("%.17g", c(n, a, prior)), collapse = "|")
  if (!is.null(hit <- .oc_cache[[key]])) return(hit)
  x <- 0:n
  A <- prior[1] + x; B <- prior[2] + n - x
  # .pr_superior_grid(): rows treatment states x1, columns control states x0
  Tab <- t(.pr_superior_grid(A, B, A, B, a, cache = FALSE))
  assign(key, Tab, envir = .oc_cache)
  Tab
}

#' Simulate a two-arm binary trial with interim looks
#'
#' Equal allocation, looks at cumulative per-arm sample sizes `looks`. Success is declared when
#' `Pr(p1 - p0 > a | data) > c` at any look; optionally the trial stops for futility at the first
#' look where that probability falls below `fut`. Efficacy stopping is applied later by
#' [oc_summarize()], so one simulation serves every candidate `c` (common random numbers) and
#' [oc_calibrate()] can find `c` from a null simulation.
#'
#' @param p0,p1 True control and treatment response rates.
#' @param looks Increasing integer vector of cumulative per-arm sample sizes.
#' @param a Superiority margin.
#' @param prior Beta prior `c(shape1, shape2)` on each arm.
#' @param nsim Number of simulated trials.
#' @param seed Optional seed.
#' @param fut Optional futility threshold on the posterior probability (binding).
#' @param method `"table"` (default, exact posterior table per look, see [oc_post_table()]) or `"normal"`
#'   (uses `prob_superior(method = "normal")`; no tables, for large `looks`).
#' @return List with `pp` (nsim x K matrix, `NA` after a futility stop), `x0`, `x1`
#'   (cumulative counts), `alive` (logical matrix, trial reaches look k), `looks`, `p0`, `p1`,
#'   `a`, `prior`.
#' @export
#' @examples
#' s <- oc_simulate(0.3, 0.3, c(50, 100), nsim = 500, seed = 1); dim(s$pp)
oc_simulate <- function(p0, p1, looks, a = 0, prior = c(1, 1), nsim = 20000, seed = NULL,
                        fut = NULL, method = c("table", "normal")) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)
  if (anyNA(looks) || any(!is.finite(looks)) || any(looks != round(looks)) || any(looks < 1))
    stop("`looks` must be positive integers")
  looks <- as.integer(looks)
  K <- length(looks)
  inc <- diff(c(0L, looks))
  draw <- function(p) {
    m <- matrix(stats::rbinom(nsim * K, rep(inc, each = nsim), p), nsim, K)
    for (k in seq_len(K)[-1]) m[, k] <- m[, k] + m[, k - 1]
    m
  }
  x0 <- draw(p0)
  x1 <- draw(p1)
  pp <- matrix(0, nsim, K)
  for (k in seq_len(K)) {
    n <- looks[k]
    if (method == "table") {
      pp[, k] <- oc_post_table(n, a, prior)[cbind(x0[, k] + 1L, x1[, k] + 1L)]
    } else {
      pp[, k] <- prob_superior(prior[1] + x1[, k], prior[2] + n - x1[, k],
                               prior[1] + x0[, k], prior[2] + n - x0[, k],
                               margin = a, method = "normal")
    }
  }
  sim <- list(pp = pp, x0 = x0, x1 = x1, alive = matrix(TRUE, nsim, K), looks = looks,
              p0 = p0, p1 = p1, a = a, prior = prior)
  if (!is.null(fut)) sim <- oc_with_futility(sim, pp, fut)
  sim
}

#' Apply a binding futility rule to a simulation
#'
#' Stops each trial at the first look where `stat < thr`; later looks are marked not alive and
#' their posterior probabilities set to `NA`.
#'
#' @param sim Output of [oc_simulate()].
#' @param stat Matrix (nsim x K) of the monitored statistic.
#' @param thr Futility threshold.
#' @return The updated `sim` list.
#' @export
#' @examples
#' s <- oc_simulate(0.3, 0.3, c(50, 100), nsim = 100, seed = 1)
#' mean(oc_with_futility(s, s$pp, 0.1)$alive[, 2])
oc_with_futility <- function(sim, stat, thr) {
  stopped <- t(apply(stat < thr, 1, cumsum)) > 0
  if (ncol(stat) == 1) stopped <- matrix(stopped, ncol = 1)
  alive <- matrix(TRUE, nrow(stat), ncol(stat))
  if (ncol(stat) > 1) alive[, -1] <- !stopped[, -ncol(stat), drop = FALSE]
  sim$alive <- alive
  sim$pp <- ifelse(alive, sim$pp, NA_real_)
  sim
}

#' Predictive probability of final success, for every (x0, x1) at a look
#'
#' `P[x0 + 1, x1 + 1]` is the probability that the rule `Pr(p1 - p0 > a | final data) > c`
#' succeeds at `Nmax` per arm, given `x0`, `x1` successes among the first `n` per arm
#' (exact beta-binomial predictive).
#'
#' @param n Per-arm patients so far.
#' @param Nmax Final per-arm sample size.
#' @param c Success threshold on the posterior probability.
#' @param a Margin.
#' @param prior Beta prior `c(shape1, shape2)` on each arm.
#' @return Matrix `(n + 1) x (n + 1)`.
#' @export
#' @examples
#' oc_predictive_table(5, 10, 0.9)[3, 4]
oc_predictive_table <- function(n, Nmax, c, a = 0, prior = c(1, 1)) {
  m <- Nmax - n
  I <- (oc_post_table(Nmax, a, prior) > c) * 1
  y <- 0:m
  x <- 0:n
  W <- t(vapply(x, function(xi) .dbetabinom(y, m, prior[1] + xi, prior[2] + n - xi),
                numeric(m + 1)))
  if (m == 0) W <- matrix(W, n + 1, 1)
  J <- t(vapply(seq_len(n + 1), function(i) as.vector(W[i, ] %*% I[i:(i + m), , drop = FALSE]),
                numeric(Nmax + 1)))
  P <- matrix(0, n + 1, n + 1)
  for (j in seq_len(n + 1)) P[, j] <- J[, j:(j + m), drop = FALSE] %*% W[j, ]
  P
}

#' Predictive probability of eventual success at every look, for every simulated trial
#'
#' @param sim Output of [oc_simulate()].
#' @param Nmax Final per-arm sample size.
#' @param c Success threshold.
#' @return Matrix (nsim x K).
#' @export
#' @examples
#' s <- oc_simulate(0.3, 0.3, c(5, 10), nsim = 50, seed = 1, method = "table")
#' dim(oc_predictive_matrix(s, 10, 0.9))
oc_predictive_matrix <- function(sim, Nmax, c) {
  out <- matrix(0, nrow(sim$pp), length(sim$looks))
  for (k in seq_along(sim$looks)) {
    tab <- oc_predictive_table(sim$looks[k], Nmax, c, sim$a, sim$prior)
    out[, k] <- tab[cbind(sim$x0[, k] + 1L, sim$x1[, k] + 1L)]
  }
  out
}

#' Largest posterior probability reached at an attended look
#'
#' @param sim Output of [oc_simulate()].
#' @return Numeric vector, one value per trial.
#' @export
#' @examples
#' s <- oc_simulate(0.3, 0.3, c(50, 100), nsim = 100, seed = 1); length(oc_max_pp(s))
oc_max_pp <- function(sim) {
  v <- ifelse(sim$alive & !is.na(sim$pp), sim$pp, -1)
  v[cbind(seq_len(nrow(v)), max.col(v, ties.method = "first"))]
}

#' Calibrate the success threshold c to a one-sided type I error
#'
#' Returns a threshold `c` at which the *empirical* type I error in the supplied null simulation
#' (the fraction of simulated trials whose posterior probability exceeds `c` at any attended
#' look) is at most `alpha`: the `1 - alpha` quantile ("higher" rule) of [oc_max_pp()]. Because
#' the rule is `pp > c`, the empirical error is at most `alpha`, with equality up to Monte Carlo
#' ties.
#'
#' The guarantee is narrower than it may look.
#' \itemize{
#'   \item It holds for the empirical rejection rate in the one null scenario that was
#'   simulated (the control rate `p0`, the looks, the margin, the prior and any futility rule in
#'   `sim_null`). It is not uniform over nuisance parameters: the same `c` can have a different
#'   type I error at another control rate, and the calibrated `c` itself varies with `p0`.
#'   \item It is a Monte Carlo estimate. A different seed or a different `nsim` gives a different
#'   `c`, and the true type I error at the returned `c` differs from `alpha` by Monte Carlo error
#'   (roughly `sqrt(alpha * (1 - alpha) / nsim)` for the empirical rate, plus the sampling
#'   variability of the quantile itself).
#'   \item The quantile convention can be one simulated trial more conservative than the
#'   smallest admissible cutoff, and `c` is a simulated value of the maximum posterior
#'   probability, so the empirical error is typically slightly below `alpha`.
#' }
#' Validate the returned `c` on fresh null simulations (new seeds, larger `nsim`) and on a
#' prespecified grid of control rates (and, if relevant, other nuisance scenarios), for example
#' with [oc_summarize()], and report the worst case or recalibrate on the least favourable rate.
#'
#' @param sim_null Output of [oc_simulate()] with `p1 = p0` (the null).
#' @param alpha One-sided type I error (default 0.025).
#' @return The calibrated threshold `c` (a number).
#' @export
#' @references Berry, S. M., Carlin, B. P., Lee, J. J., Muller, P. (2010). Bayesian Adaptive
#'   Methods for Clinical Trials. CRC Press.
#' @examples
#' oc_calibrate(oc_simulate(0.3, 0.3, c(50, 100), nsim = 2000, seed = 1))
oc_calibrate <- function(sim_null, alpha = 0.025) {
  m <- sort(oc_max_pp(sim_null))
  n <- length(m)
  idx <- ceiling((n - 1) * (1 - alpha) - 1e-9) + 1
  m[min(idx, n)]
}

#' Operating characteristics of the rule "success if pp > c at any attended look"
#'
#' @param sim Output of [oc_simulate()].
#' @param c Success threshold.
#' @return Named list: `reject` (power or type I error), `EN` (expected total sample size,
#'   both arms), `early_succ`, `fut_stop`, `bias` (all trials), `bias_succ`, `bias_early`
#'   (bias of the observed difference among successes and early successes) and `se_reject`
#'   (Monte Carlo standard error of `reject`).
#' @export
#' @examples
#' s <- oc_simulate(0.3, 0.5, c(50, 100), nsim = 500, seed = 1); oc_summarize(s, 0.95)$reject
oc_summarize <- function(sim, c) {
  pp <- sim$pp; alive <- sim$alive; looks <- sim$looks
  nsim <- nrow(pp); K <- ncol(pp)
  win <- (pp > c) & alive
  win[is.na(win)] <- FALSE
  anyw <- rowSums(win) > 0
  first <- ifelse(anyw, max.col(win, ties.method = "first"), K + 1L)
  success <- first <= K
  att <- rowSums(alive)
  stop_k <- ifelse(success, first, att)
  N <- 2 * looks[stop_k]
  early <- success & (first < K)
  fut_stop <- !success & (att < K)
  idx <- cbind(seq_len(nsim), stop_k)
  n <- looks[stop_k]
  est <- sim$x1[idx] / n - sim$x0[idx] / n
  truth <- sim$p1 - sim$p0
  ps <- mean(success)
  list(reject = ps, EN = mean(N), early_succ = mean(early), fut_stop = mean(fut_stop),
       bias = mean(est - truth),
       bias_succ = if (any(success)) mean(est[success] - truth) else NA_real_,
       bias_early = if (any(early)) mean(est[early] - truth) else NA_real_,
       se_reject = sqrt(ps * (1 - ps) / nsim))
}

.dbetabinom <- function(x, m, a, b) {
  exp(lchoose(m, x) + lbeta(x + a, m - x + b) - lbeta(a, b))
}
