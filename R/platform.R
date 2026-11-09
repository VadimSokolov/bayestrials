# Master protocols and platform trials: helper calculations and simulations.

#' Shared-control enrolment table
#'
#' Total enrolment for separate trials versus one platform with a shared control of size `m`.
#'
#' @param ks Numbers of experimental arms.
#' @param m Patients per arm.
#' @return Data frame with `k`, `separate`, `platform`, `saving`.
#' @examples
#' platform_shared_control_table(2:6)
#' @export
platform_shared_control_table <- function(ks, m = 80) {
  sep <- 2 * ks * m; plat <- (ks + 1) * m
  data.frame(k = ks, separate = sep, platform = plat, saving = 1 - plat / sep)
}

#' Variance of one arm-vs-control mean difference
#'
#' Sigma = 1, `N` patients split with `ratio` control patients per patient on each of `k` arms.
#'
#' @param k Number of experimental arms.
#' @param N Total patients.
#' @param ratio Control to arm allocation ratio.
#' @return Variance.
#' @examples
#' platform_contrast_variance(4, 480, 2)
#' @export
platform_contrast_variance <- function(k, N, ratio) {
  n <- N / (ratio + k)
  1 / n + 1 / (ratio * n)
}

#' Correlation between two arm-vs-control z statistics
#'
#' @param ratio Control to arm allocation ratio.
#' @return Correlation `1 / (1 + ratio)`.
#' @examples
#' platform_arm_correlation(2)
#' @export
platform_arm_correlation <- function(ratio) 1 / (1 + ratio)

#' Distribution of the number of false positives with a shared control
#'
#' All arms null, each of `k` arm-vs-control one-sided tests at level `alpha`.
#'
#' @param k Number of arms.
#' @param ratio Control to arm allocation ratio.
#' @param alpha One-sided level.
#' @param nsim Number of simulated platforms.
#' @param seed Optional seed.
#' @return Numeric vector of length `k + 1`: probabilities of 0, ..., k false positives.
#' @examples
#' platform_false_positive_counts(3, 1, nsim = 1000, seed = 1)
#' @export
platform_false_positive_counts <- function(k, ratio, alpha = 0.025, nsim = 400000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- 100; n0 <- ratio * n
  cc <- rnorm(nsim, 0, 1 / sqrt(n0))
  a <- matrix(rnorm(nsim * k, 0, 1 / sqrt(n)), nsim, k)
  z <- (a - cc) / sqrt(1 / n + 1 / n0)
  fp <- rowSums(z > qnorm(alpha, lower.tail = FALSE))
  tabulate(fp + 1, nbins = k + 1) / nsim
}

#' Non-concurrent control drift simulation
#'
#' Arm A opens at t = 0 and arm B at `b_open` years. Control success rises by `drift` per year.
#' Arm B is analysed by concurrent controls, all controls, and a time-adjusted regression.
#'
#' @param drift Control drift per year.
#' @param effect True effect of arm B.
#' @param p0 Control rate at time 0.
#' @param per_year Patients per year.
#' @param years Platform duration.
#' @param b_open Opening time of arm B.
#' @param nsim Number of simulations.
#' @param seed Optional seed.
#' @return Data frame with `method`, `mean_estimate`, `rejection`.
#' @examples
#' platform_drift_simulation(nsim = 200, seed = 1)
#' @export
platform_drift_simulation <- function(drift = 0.02, effect = 0, p0 = 0.30, per_year = 300,
                                      years = 4, b_open = 2, nsim = 20000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  zc <- qnorm(0.025, lower.tail = FALSE)
  t_all <- numeric(0); g_all <- integer(0)
  for (yr in 0:(years - 1)) {
    for (i in 0:(per_year - 1)) {
      t <- yr + (i + 0.5) / per_year
      arms <- if (t < b_open) c(0, 1) else c(0, 1, 2)
      g_all <- c(g_all, arms[i %% length(arms) + 1])
      t_all <- c(t_all, t)
    }
  }
  ctrl <- g_all == 0; isB <- g_all == 2; conc <- ctrl & t_all >= b_open
  pc <- p0 + drift * t_all
  p <- ifelse(isB, pc + effect, pc)
  N <- length(p)
  sel <- ctrl | isB
  X <- cbind(1, as.numeric(isB[sel]), t_all[sel])
  XtXi <- solve(crossprod(X)); H <- XtXi %*% t(X)
  nb <- sum(isB); nct <- sum(conc); nal <- sum(ctrl)
  est <- matrix(0, nsim, 3); rej <- matrix(FALSE, nsim, 3)
  chunk <- 2000
  start <- 1
  while (start <= nsim) {
    m <- min(chunk, nsim - start + 1)
    y <- matrix(as.numeric(runif(m * N) < rep(p, each = m)), m, N)
    yb <- y[, isB, drop = FALSE]
    mb <- rowMeans(yb); vb <- rowSums((yb - mb)^2) / (nb - 1)
    for (s in 1:2) {
      yc <- if (s == 1) y[, conc, drop = FALSE] else y[, ctrl, drop = FALSE]
      nn <- if (s == 1) nct else nal
      mc <- rowMeans(yc); vc <- rowSums((yc - mc)^2) / (nn - 1)
      e <- mb - mc; se <- sqrt(vb / nb + vc / nn)
      est[start:(start + m - 1), s] <- e
      rej[start:(start + m - 1), s] <- e / se > zc
    }
    yy <- y[, sel, drop = FALSE]
    beta <- yy %*% t(H)
    rss <- rowSums(yy^2) - rowSums(beta * (yy %*% X))
    se <- sqrt(XtXi[2, 2] * rss / (ncol(yy) - 3))
    est[start:(start + m - 1), 3] <- beta[, 2]
    rej[start:(start + m - 1), 3] <- beta[, 2] / se > zc
    start <- start + m
  }
  data.frame(method = c("concurrent", "all controls", "time-adjusted"),
             mean_estimate = colMeans(est), rejection = colMeans(rej))
}

#' Predictive probability of phase 3 success (graduation)
#'
#' Phase 3 with `n3` per arm, 1:1, one-sided 0.025 z test, given phase 2 counts, Beta(1, 1) priors.
#'
#' @param x_t,n_t Phase 2 treatment responders and patients.
#' @param x_c,n_c Phase 2 control responders and patients.
#' @param n3 Phase 3 patients per arm.
#' @param nsim Number of simulations.
#' @param seed Optional seed.
#' @return Named vector: `pp` (predictive probability of phase 3 success) and `post` (posterior probability treatment is better).
#' @examples
#' platform_predictive_graduation(20, 40, 12, 40, nsim = 2000, seed = 1)
#' @export
platform_predictive_graduation <- function(x_t, n_t, x_c, n_c, n3 = 150, nsim = 40000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  pt <- rbeta(nsim, 1 + x_t, 1 + n_t - x_t)
  pcn <- rbeta(nsim, 1 + x_c, 1 + n_c - x_c)
  xt3 <- rbinom(nsim, n3, pt); xc3 <- rbinom(nsim, n3, pcn)
  pool <- (xt3 + xc3) / (2 * n3)
  se <- sqrt(pmax(pool * (1 - pool) * 2 / n3, 1e-12))
  z <- (xt3 / n3 - xc3 / n3) / se
  c(pp = mean(z > qnorm(0.025, lower.tail = FALSE)), post = mean(pt > pcn))
}

# Smoothing fit for one data set: uses only the observed counts `y` (rows = replicates, columns =
# periods) and the per-period sample size `n_per`. Nothing about the data-generating rates enters.
# Observation variances are estimated from the data, v_j = pt_j (1 - pt_j) / n_per with
# pt_j = (y_j + 0.5) / (n_per + 1). The period means follow a first-order Gaussian random walk,
# mu_j = mu_{j-1} + N(0, tau^2), with a flat prior on mu_1. Given tau the estimate of the last period
# mean is the posterior mean. With `taus = NULL`, tau is estimated per replicate by maximising the
# marginal (restricted) likelihood of the period differences, Dy ~ N(0, D V D' + tau^2 I).
.tm_fit <- function(y, n_per, taus = NULL, tau_grid = exp(seq(log(1e-3), log(1), length.out = 40))) {
  y <- as.matrix(y)
  J <- ncol(y)
  ybar <- y / n_per
  if (J == 1L) {
    # a single period: nothing to smooth, tau is not identified
    if (!is.null(taus)) {
      out <- lapply(taus, function(tau) ybar[, 1]); names(out) <- paste0("fixed smoothing, tau=", taus)
      return(list(est = out, tau_hat = NULL))
    }
    return(list(est = list(`learned smoothing (REML tau)` = ybar[, 1]), tau_hat = rep(NA_real_, nrow(y))))
  }
  pt <- (y + 0.5) / (n_per + 1)
  v <- pt * (1 - pt) / n_per
  D <- diff(diag(J))
  DtD <- crossprod(D)
  wts_for <- function(vv, tau) {
    P <- diag(1 / vv, J) + DtD / tau^2
    solve(P)[J, ] / vv
  }
  if (!is.null(taus)) {
    out <- lapply(taus, function(tau) {
      vapply(seq_len(nrow(y)), function(r) sum(wts_for(v[r, ], tau) * ybar[r, ]), numeric(1))
    })
    names(out) <- paste0("fixed smoothing, tau=", taus)
    return(list(est = out, tau_hat = NULL))
  }
  dy <- ybar %*% t(D)
  tau_hat <- numeric(nrow(y)); est <- numeric(nrow(y))
  for (r in seq_len(nrow(y))) {
    DVD <- D %*% (v[r, ] * t(D))
    ll <- vapply(tau_grid, function(tau) {
      Q <- DVD + diag(tau^2, J - 1)
      ch <- chol(Q)
      -sum(log(diag(ch))) - 0.5 * sum(backsolve(ch, dy[r, ], transpose = TRUE)^2)
    }, numeric(1))
    tau_hat[r] <- tau_grid[which.max(ll)]
    est[r] <- sum(wts_for(v[r, ], tau_hat[r]) * ybar[r, ])
  }
  list(est = list(`learned smoothing (REML tau)` = est), tau_hat = tau_hat)
}

#' Time machine: estimating the current control rate under period drift
#'
#' Simulation study of estimators of the control rate in the current (last) period when the true control rate drifts across periods. The true rate in period `j` is `p0 + shift[j]` and is used only to generate the binomial counts and to score the estimators. Every estimator is computed from the simulated counts and `n_per` alone, so identical observed counts give identical estimates whatever the scenario that produced them.
#'
#' Estimators: `current only` (last period mean), `pool all` (mean over all periods), `learned smoothing (REML tau)`, and `fixed smoothing, tau=...` for each value in `taus`. The smoothers are first-order Gaussian random walks on the period means with step sd `tau`, a flat prior on the first period, and observation variances estimated from the observed period proportions (with a 0.5 continuity correction). In the learned version `tau` is estimated in each replicate by maximising the marginal (restricted) likelihood of the period differences over a log-spaced grid on `[0.001, 1]`; the fixed versions are illustrations with a prespecified `tau`. This is a first-order Gaussian approximation in the spirit of the Bayesian Time Machine, not the published second-order binomial model: it does not reproduce its results.
#'
#' @param shift Vector of period shifts (length `periods`).
#' @param nsim Number of simulations.
#' @param seed Optional seed.
#' @param p0 Baseline rate.
#' @param n_per Patients per period.
#' @param taus Random-walk step sds for the fixed-smoothing illustrations.
#' @return Data frame with `method`, `bias`, `rmse`, and `mean_tau` (average step sd used: the average estimate for the learned smoother, `tau` for fixed smoothers, `NA` otherwise).
#' @examples
#' platform_time_machine(c(0, 0.02, 0.04, 0.06), nsim = 200, seed = 1)
#' @export
platform_time_machine <- function(shift, nsim = 4000, seed = NULL, p0 = 0.30, n_per = 100,
                                  taus = c(0.01, 0.03, 0.10)) {
  if (!is.null(seed)) set.seed(seed)
  truth <- p0 + as.numeric(shift)
  J <- length(truth)
  y <- matrix(rbinom(nsim * J, n_per, rep(truth, each = nsim)), nsim, J)
  # from here on only the observed counts are used
  ybar <- y / n_per
  learned <- .tm_fit(y, n_per)
  fixed <- .tm_fit(y, n_per, taus = taus)
  est <- c(list(`current only` = ybar[, J], `pool all` = rowMeans(ybar)), learned$est, fixed$est)
  mean_tau <- c(NA, NA, if (J == 1L) NA_real_ else mean(learned$tau_hat), taus)
  tj <- truth[J]
  data.frame(method = names(est),
             bias = vapply(est, function(e) mean(e) - tj, numeric(1)),
             rmse = vapply(est, function(e) sqrt(mean((e - tj)^2)), numeric(1)),
             mean_tau = mean_tau,
             row.names = NULL)
}
