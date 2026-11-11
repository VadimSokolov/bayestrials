# Decision-theoretic design tools for the book: Bayes risk of a test, Jeffreys bar in q,
# square-root sizing, cost of K subsets, Bayesian sample size criteria, optimal stopping,
# pivotal sizing, noninferiority and borrowing worked cases.
# Notation: q = I tau^2 is the prior-to-sampling information ratio of the normal mixture
# (H1: theta ~ N(0, tau^2), Z | H1 ~ N(0, 1 + q)); the Bayes factor is
# B = (1 + q)^(-1/2) exp(z^2 q / (2 (1 + q))).

#' Bayes factor of the normal mixture against the point null, in terms of q
#'
#' \eqn{B = (1+q)^{-1/2}\exp\{z^2 q/(2(1+q))\}} for \eqn{H_1:\theta\sim N(0,\tau^2)},
#' \eqn{q = I\tau^2}.
#'
#' @param z Standardized statistic.
#' @param q Information ratio \eqn{q = I\tau^2} (for a unit-information prior, the sample size).
#' @return Numeric vector of Bayes factors in favor of the effect.
#' @examples
#' decide_bf_q(1.96, 10)
#' @export
decide_bf_q <- function(z, q) exp(-0.5 * log1p(q) + z^2 * q / (2 * (1 + q)))

#' Exact Bayes evidence bar for the normal mixture
#'
#' The statistic bar at which \eqn{B = c} is
#' \eqn{b = \sqrt{(1+1/q)(\log(1+q) + 2\log c)}}. This is the exact bar; the large-information
#' form \eqn{\sqrt{\log(1+q)}} drops the factor \eqn{(1+1/q)} and is not exact.
#'
#' @param q Information ratio \eqn{q = I\tau^2}.
#' @param c Bayes factor threshold (default 1, i.e. equal prior odds and equal losses).
#' @return Numeric vector of bars on the \eqn{|Z|} scale. A threshold below the minimum attainable Bayes factor \eqn{(1+q)^{-1/2}} makes rejection certain, and the bar is 0.
#' @examples
#' decide_bayes_bar(c(10, 100))
#' @export
decide_bayes_bar <- function(q, c = 1) sqrt((1 + 1 / q) * pmax(log1p(q) + 2 * log(c), 0))

#' Information ratio at which the exact Bayes bar equals a given bar
#'
#' Solves \eqn{(1+1/q)\log(1+q) = z^2}; for \eqn{z = 1.96} the root is \eqn{q\approx 41.6}.
#' A positive root exists only when \eqn{|z| > 1} (the left side increases from 1 at \eqn{q=0}).
#'
#' @param z Bar (default 1.96).
#' @return The crossing value of \eqn{q}.
#' @examples
#' decide_crossing_q()
#' @export
decide_crossing_q <- function(z = 1.96) {
  z <- abs(z)
  if (!(z > 1)) stop("a positive crossing exists only for abs(z) > 1")
  f <- function(q) (1 + 1 / q) * log1p(q) - z^2
  lo <- 1e-10; hi <- 2
  while (f(hi) < 0) { hi <- hi * 2; if (hi > 1e300) stop("no crossing found") }
  stats::uniroot(f, c(lo, hi), tol = 1e-10)$root
}

#' Bayes risk of a two-sided test in the normal mixture model
#'
#' The test rejects when \eqn{|Z| > b}. Under \eqn{H_0} \eqn{Z\sim N(0,1)}; under \eqn{H_1}
#' \eqn{Z\sim N(0, 1+q)}. The Bayes risk is
#' \eqn{\pi_0\ell_1\alpha + (1-\pi_0)\ell_0\beta}.
#'
#' @param b Rejection bar on the \eqn{|Z|} scale (default 1.96, the fixed-level test).
#' @param q Information ratio \eqn{q = I\tau^2}.
#' @param pi0 Prior probability of the null.
#' @param l1,l0 Loss of a false approval and of a false rejection.
#' @return Data frame with columns `b`, `q`, `alpha`, `beta`, `risk`.
#' @examples
#' decide_bayes_risk(1.96, c(10, 1000))
#' @export
decide_bayes_risk <- function(b = 1.96, q, pi0 = 0.5, l1 = 1, l0 = 1) {
  n <- max(length(b), length(q)); b <- rep_len(b, n); q <- rep_len(q, n)
  alpha <- 2 * stats::pnorm(-b)
  beta <- 2 * stats::pnorm(b / sqrt(1 + q)) - 1
  data.frame(b = b, q = q, alpha = alpha, beta = beta,
             risk = pi0 * l1 * alpha + (1 - pi0) * l0 * beta)
}

#' Bayes factor threshold implied by a fixed bar
#'
#' The test \eqn{|Z|>b} is the Bayes rule \eqn{B>c} with
#' \eqn{2\log c = b^2 q/(1+q) - \log(1+q)}, exactly. The product of prior odds and loss ratio
#' that makes it optimal is \eqn{c}.
#'
#' @param b Bar on the \eqn{|Z|} scale.
#' @param q Information ratio.
#' @return Numeric vector of implied thresholds \eqn{c}.
#' @examples
#' decide_implied_c(1.96, c(10, 1000))
#' @export
decide_implied_c <- function(b, q) exp((b^2 * q / (1 + q) - log(1 + q)) / 2)

#' Loss ratio implied by a posterior probability threshold, and its inverse
#'
#' Acting when \eqn{\Pr(H_1\mid y) > p^\star} is optimal when a false approval is
#' \eqn{p^\star/(1-p^\star)} times as costly as a false rejection (equal prior odds are not
#' needed: the ratio is of posterior odds).
#'
#' @param pstar Posterior probability threshold, or NULL.
#' @param ratio Loss ratio \eqn{\ell_1/\ell_0}, or NULL.
#' @return Data frame with `pstar` and `ratio`.
#' @examples
#' decide_loss_ratio(c(0.975, 0.99, 0.995))
#' @export
decide_loss_ratio <- function(pstar = NULL, ratio = NULL) {
  if (!is.null(pstar)) ratio <- pstar / (1 - pstar) else pstar <- ratio / (1 + ratio)
  data.frame(pstar = pstar, ratio = ratio)
}

#' Lindley paradox table
#'
#' Bayes factor and posterior probability of the null at a fixed statistic \eqn{Z} as the
#' trial grows, with prior \eqn{\tau} on the scale of one observation.
#'
#' @param n Sample sizes.
#' @param z Observed statistic (default 1.96).
#' @param tau Prior standard deviation of the effect under \eqn{H_1}.
#' @param sigma Observation standard deviation.
#' @param pi0 Prior probability of the null.
#' @return Data frame with `n`, `q`, `B` and `post_null`.
#' @examples
#' decide_lindley_table()
#' @export
decide_lindley_table <- function(n = c(10, 100, 1000, 1e4, 1e5), z = 1.96, tau = 1, sigma = 1,
                                 pi0 = 0.5) {
  q <- n * tau^2 / sigma^2
  B <- decide_bf_q(z, q)
  data.frame(n = n, q = q, B = B, post_null = pi0 / (pi0 + (1 - pi0) * B))
}

#' Type I error of the Bayes rule
#'
#' Two-sided type I error of the rule \eqn{B>c}, reported with the exact bar
#' \eqn{\sqrt{(1+1/q)\log(1+q)}} and, for comparison, the large-information bar
#' \eqn{\sqrt{\log(1+q)}} that the zip book prints.
#'
#' @param n Sample sizes (with a unit information prior, \eqn{q = n\tau^2/\sigma^2}).
#' @param tau,sigma Prior sd and observation sd.
#' @param c Bayes factor threshold.
#' @return Data frame with `n`, `q`, `bar_exact`, `alpha_exact`, `bar_approx`, `alpha_approx`.
#' @examples
#' decide_bayes_type1_table()
#' @export
decide_bayes_type1_table <- function(n = c(10, 100, 1000, 1e4), tau = 1, sigma = 1, c = 1) {
  q <- n * tau^2 / sigma^2
  be <- decide_bayes_bar(q, c); ba <- sqrt(pmax(log1p(q) + 2 * log(c), 0))
  data.frame(n = n, q = q, bar_exact = be, alpha_exact = 2 * stats::pnorm(-be),
             bar_approx = ba, alpha_approx = 2 * stats::pnorm(-ba))
}

#' Exact Bayes factor in favor of the null under Jeffreys' Cauchy prior
#'
#' Prior \eqn{\delta\sim}Cauchy(0,1) on the effect in units of \eqn{\sigma}; the standardized
#' effect \eqn{\sqrt n\delta} has scale \eqn{\sqrt n}.
#'
#' @param t Standardized statistic.
#' @param n Sample size.
#' @return Bayes factor \eqn{K = p(t\mid H_0)/p(t\mid H_1)}.
#' @examples
#' decide_jeffreys_k_exact(sqrt(log(pi * 10 / 2)), 10)
#' @export
decide_jeffreys_k_exact <- function(t, n) {
  len <- max(length(t), length(n)); t <- rep_len(t, len); n <- rep_len(n, len)
  vapply(seq_len(len), function(i) {
    # centered substitution u = l - t puts the normal peak at zero
    m1 <- stats::integrate(function(u) stats::dnorm(u) * stats::dcauchy(t[i] + u, 0, sqrt(n[i])),
                           -Inf, Inf, subdivisions = 2000L, rel.tol = 1e-10)$value
    exp(stats::dnorm(t[i], log = TRUE) - log(m1))
  }, numeric(1))
}

#' Jeffreys' evidence table, with the normal-mixture analogue in q
#'
#' Columns for the Cauchy large-n criterion \eqn{t=\sqrt{\log(\pi n/2)}}, its two sided p value,
#' the exact Cauchy Bayes factor there, the \eqn{t} for \eqn{K^{-1}=10}, and the exact
#' normal-mixture bars (Bayes factor 1 and 10) with \eqn{q = n\tau^2}, which carry the factor
#' \eqn{(1+1/q)}.
#'
#' @param n Sample sizes.
#' @param tau Prior sd (unit observation sd) for the normal-mixture columns.
#' @param strong Bayes factor against the null for the strong-evidence column. `t_strong` is the large-sample
#'   approximation \eqn{\sqrt{\log(K^2\pi n/2)}}, not the exact solution of the Bayes factor equation (at \eqn{n=10} it gives an
#'   exact Bayes factor of about 5.9, not 10).
#' @return Data frame.
#' @examples
#' decide_jeffreys_table()
#' @export
decide_jeffreys_table <- function(n = c(10, 100, 1e3, 1e4, 3e4, 1e5, 1e6), tau = 1, strong = 10) {
  tc <- sqrt(log(pi * n / 2)); q <- n * tau^2
  data.frame(n = n, t_crit = tc, p_two_sided = 2 * stats::pnorm(-tc),
             K_exact_at_tcrit = decide_jeffreys_k_exact(tc, n),
             t_strong = sqrt(log(strong^2 * pi * n / 2)),
             q = q, bar_mix_1 = decide_bayes_bar(q, 1),
             bar_mix_strong = decide_bayes_bar(q, strong), conventional = 1.96)
}

# ---- square-root sizing --------------------------------------------------------------------

#' Total expected loss of a two-arm trial of n per arm (planner's problem)
#'
#' \eqn{L(n) = n\tau\sqrt{2/\pi} + (N-2n)\tau(1-\rho)/\sqrt{2\pi} + 2cn}, with
#' \eqn{\rho = \tau/\sqrt{\tau^2+2\sigma^2/n}}.
#'
#' @param n Patients per arm.
#' @param N Population size.
#' @param tau Prior sd of the effect (prior mean 0).
#' @param sigma Outcome sd.
#' @param cost Monetary cost per enrolled patient, in loss units.
#' @return Numeric vector of losses.
#' @examples
#' size_sqrt_loss(1084, 1e6, 1.3, 2)
#' @export
size_sqrt_loss <- function(n, N, tau, sigma, cost = 0) {
  rho <- tau / sqrt(tau^2 + 2 * sigma^2 / n)
  n * tau * sqrt(2 / pi) + (N - 2 * n) * tau * (1 - rho) / sqrt(2 * pi) + 2 * cost * n
}

#' Square-root sizing rule
#'
#' \eqn{n^\star = \sigma\sqrt{N/(2\tau^2 + 2\sqrt{2\pi}c\tau)}} per arm (\eqn{c=0}: the
#' Anscombe-Colton rule \eqn{(\sigma/\tau)\sqrt{N/2}}), with approximate minimized loss
#' \eqn{2\sigma\sqrt{N/\pi}} at \eqn{c=0}.
#'
#' @inheritParams size_sqrt_loss
#' @return Named numeric vector `c(n_star, loss_approx)`.
#' @examples
#' size_sqrt_rule(1e6, 1.3, 2)
#' @export
size_sqrt_rule <- function(N, tau, sigma, cost = 0) {
  n <- sigma * sqrt(N / (2 * tau^2 + 2 * sqrt(2 * pi) * cost * tau))
  c(n_star = n, loss_approx = if (cost == 0) 2 * sigma * sqrt(N / pi) else NA_real_)
}

#' Numerical optimizer of the square-root loss, compared with the closed form
#'
#' @param N Vector of population sizes.
#' @inheritParams size_sqrt_loss
#' @return Data frame with `N`, `n_star`, `n_opt`, `loss_approx`, `loss_min`, `rel_error`.
#' @examples
#' size_sqrt_table(c(1e5, 1e6), 1.3, 2)
#' @export
size_sqrt_table <- function(N = c(1e5, 1e6), tau = 1.3, sigma = 2, cost = 0) {
  do.call(rbind, lapply(N, function(Ni) {
    r <- size_sqrt_rule(Ni, tau, sigma, cost)
    o <- stats::optimize(function(n) size_sqrt_loss(n, Ni, tau, sigma, cost),
                         c(1, Ni / 2), tol = 1e-9)
    # integer minimizer for reporting
    cand <- floor(o$minimum) + (-1:2)
    Lc <- size_sqrt_loss(cand, Ni, tau, sigma, cost)
    data.frame(N = Ni, n_star = unname(r[1]), n_opt = cand[which.min(Lc)],
               loss_approx = unname(r[2]), loss_min = min(Lc),
               rel_error = if (cost == 0) (unname(r[2]) - min(Lc)) / min(Lc) else NA_real_)
  }))
}

# ---- cost of K subsets and the max-of-K simulation -----------------------------------------

#' Cost of preserving K analyses (Sidak threshold and extra patients)
#'
#' @param K Numbers of analyses.
#' @param alpha One-sided level.
#' @param power Target power.
#' @return Data frame with `K`, `z_K`, `ratio` (\eqn{n_K/n_1}), `extra_per_1000`, `logK`.
#' @examples
#' decide_subset_cost()
#' @export
decide_subset_cost <- function(K = c(1, 5, 10, 50, 100, 1000), alpha = 0.025, power = 0.8) {
  zK <- stats::qnorm((1 - alpha)^(1 / K)); zb <- stats::qnorm(power)
  r <- ((zK + zb) / (stats::qnorm(1 - alpha) + zb))^2
  data.frame(K = K, z_K = zK, ratio = r, extra_per_1000 = 1000 * (r - 1), logK = log(K))
}

#' Selected-likelihood max-of-K simulation
#'
#' Compares three rules on K equicorrelated standardized statistics: the naive maximum tested
#' at one-sided \eqn{\alpha}, the Sidak maximum, and the equal-weight mixture e value
#' \eqn{\bar e = K^{-1}\sum_k (1+\lambda^2)^{-1/2}\exp\{\lambda^2 Z_k^2/(2(1+\lambda^2))\}}
#' at \eqn{\bar e\ge 1/\alpha}. Power is computed at thresholds recalibrated to realized size
#' \eqn{\alpha}. For \eqn{K=20}, \eqn{\lambda^2=9} the mixture loses power to Sidak (a
#' negative result for the mixture).
#'
#' @param K Number of statistics.
#' @param rho Vector of equicorrelations for the size panel (power uses the first two).
#' @param alpha One-sided level.
#' @param lam2 Prior variance \eqn{\lambda^2} of the mixture.
#' @param nsim Replicates (the zip uses 4e5).
#' @param scenarios Named list of effect vectors (non-null leading entries), default the zip's.
#' @param power_rho Correlations at which power is computed.
#' @param seed Seed or NULL.
#' @param direction `"two.sided"` (default) uses the symmetric squared mixture above. `"positive"` uses the
#'   positive-direction (one-sided) mixture \eqn{K^{-1}\sum_k 2(1+\lambda^2)^{-1/2}\exp\{\lambda^2 Z_k^2/(2(1+\lambda^2))\}\mathbb{1}\{Z_k>0\}},
#'   which tests the same alternative as the maximum.
#' @return List with data frames `size` (naive, sidak, mixture), `power` (sidak, mixture) and
#'   `thresholds`.
#' @examples
#' decide_max_of_k_sim(K = 5, rho = 0, nsim = 2000, scenarios = list("1 effect" = 3), power_rho = 0, seed = 1)
#' @export
decide_max_of_k_sim <- function(K = 20, rho = c(0, 0.5, 0.8), alpha = 0.025, lam2 = 9,
                                nsim = 2e5,
                                scenarios = list("one effect d=3.0" = 3, "one effect d=3.5" = 3.5,
                                                 "three effects d=2.5" = rep(2.5, 3),
                                                 "five effects d=2.2" = rep(2.2, 5)),
                                power_rho = c(0, 0.5), seed = NULL, direction = c("two.sided", "positive")) {
  direction <- match.arg(direction)
  if (!is.null(seed)) set.seed(seed)
  c0 <- (1 + lam2)^-0.5; bb <- lam2 / (2 * (1 + lam2))
  draw <- function(theta, r) {
    common <- stats::rnorm(nsim)
    z <- sqrt(r) * common + sqrt(1 - r) * matrix(stats::rnorm(nsim * K), nsim, K)
    sweep(z, 2, theta, "+")
  }
  stat <- function(z) list(mx = do.call(pmax, as.data.frame(z)),
                           E = if (direction == "positive") rowMeans(2 * c0 * exp(bb * z^2) * (z > 0))
                               else rowMeans(c0 * exp(bb * z^2)))
  zero <- numeric(K); size <- list(); cal <- list()
  for (r in rho) {
    s <- stat(draw(zero, r))
    size[[length(size) + 1]] <- data.frame(
      rho = r, naive = mean(s$mx > stats::qnorm(1 - alpha)),
      sidak = mean(s$mx > stats::qnorm((1 - alpha)^(1 / K))),
      mixture = mean(s$E >= 1 / alpha))
    cal[[as.character(r)]] <- c(sidak = unname(stats::quantile(s$mx, 1 - alpha)),
                                mixture = unname(stats::quantile(s$E, 1 - alpha)))
  }
  pw <- list()
  for (r in power_rho) {
    th <- cal[[as.character(r)]]
    if (is.null(th)) {
      s <- stat(draw(zero, r))
      th <- c(sidak = unname(stats::quantile(s$mx, 1 - alpha)),
              mixture = unname(stats::quantile(s$E, 1 - alpha)))
    }
    for (nm in names(scenarios)) {
      t <- zero; t[seq_along(scenarios[[nm]])] <- scenarios[[nm]]
      s <- stat(draw(t, r))
      pw[[length(pw) + 1]] <- data.frame(scenario = nm, rho = r,
                                         sidak = mean(s$mx > th[["sidak"]]),
                                         mixture = mean(s$E >= th[["mixture"]]))
    }
  }
  list(size = do.call(rbind, size), power = do.call(rbind, pw),
       thresholds = do.call(rbind, lapply(names(cal), function(r)
         data.frame(rho = as.numeric(r), sidak = cal[[r]][["sidak"]],
                    mixture = cal[[r]][["mixture"]]))))
}

# ---- Bayesian sample size criteria (Inoue, Berry, Parmigiani 2005) -----------------------------

# posterior summary for normal mean, two arms of n, known sigma; analysis prior N(ma, sa^2) (sa = Inf flat)
.ss_post <- function(n, sigma, ma, sa) {
  v <- 2 * sigma^2 / n
  prec <- 1 / v + if (is.finite(sa)) 1 / sa^2 else 0
  list(v = v, pv = 1 / prec, w = (1 / v) / prec, ma = if (is.finite(sa)) ma else 0, sa = sa)
}

#' Average length criterion (ALC): smallest n per arm giving a posterior interval of length l
#'
#' Normal outcomes, known variance, two arms of n. The posterior sd does not depend on the data,
#' so the expected length is \eqn{2z_{1-\alpha/2}\sqrt{\text{post. var}}}.
#'
#' @param l Target interval length.
#' @param sigma Outcome sd.
#' @param alpha Total non-coverage (default 0.05).
#' @param analysis_sd Analysis prior sd (`Inf` = flat).
#' @return Smallest integer n per arm.
#' @examples
#' decide_ss_alc(0.3)
#' @export
decide_ss_alc <- function(l, sigma = 1, alpha = 0.05, analysis_sd = Inf) {
  z <- stats::qnorm(1 - alpha / 2)
  pv_target <- (l / (2 * z))^2
  prec_needed <- 1 / pv_target - if (is.finite(analysis_sd)) 1 / analysis_sd^2 else 0
  if (prec_needed <= 0) return(0)
  ceiling(2 * sigma^2 * prec_needed)
}

# true (design-prior) coverage of the analysis interval, as function of xbar
.ss_cover <- function(xbar, n, sigma, l, dm, ds, ma, sa) {
  p <- .ss_post(n, sigma, ma, sa)
  mua <- p$w * xbar + (1 - p$w) * p$ma
  pd <- .ss_post(n, sigma, dm, ds)
  mud <- pd$w * xbar + (1 - pd$w) * pd$ma
  s <- sqrt(pd$pv)
  stats::pnorm((mua + l / 2 - mud) / s) - stats::pnorm((mua - l / 2 - mud) / s)
}

#' Average coverage (ACC) and worst outcome (WOC) criteria for a fixed-length interval
#'
#' The interval has fixed length l and is centred at the analysis-posterior mean. Coverage is
#' computed under the design prior \eqn{N(m_d, s_d^2)}. ACC averages the coverage over the
#' marginal distribution of the data; WOC is the minimum coverage over the central `gamma`
#' region of the data (the worst outcome in a set of probability `gamma`). With matched design
#' and analysis priors and known variance the coverage is constant in the data and ACC = WOC.
#'
#' @param n Patients per arm.
#' @param l Interval length.
#' @param sigma Outcome sd.
#' @param design_mean,design_sd Design prior.
#' @param analysis_mean,analysis_sd Analysis prior (`analysis_sd = Inf`: flat).
#' @param gamma Probability of the data region for WOC.
#' @return Data frame with `n`, `ACC`, `WOC`.
#' @examples
#' decide_ss_acc_woc(100, 0.4)
#' @export
decide_ss_acc_woc <- function(n, l, sigma = 1, design_mean = 0.3, design_sd = 0.1,
                              analysis_mean = 0, analysis_sd = Inf, gamma = 0.9) {
  do.call(rbind, lapply(n, function(ni) {
    sx <- sqrt(design_sd^2 + 2 * sigma^2 / ni)
    f <- function(x) .ss_cover(x, ni, sigma, l, design_mean, design_sd, analysis_mean, analysis_sd)
    acc <- stats::integrate(function(x) f(x) * stats::dnorm(x, design_mean, sx),
                            design_mean - 10 * sx, design_mean + 10 * sx)$value
    xs <- design_mean + sx * stats::qnorm(c((1 - gamma) / 2, 1 - (1 - gamma) / 2))
    xg <- seq(xs[1], xs[2], length.out = 201)
    data.frame(n = ni, ACC = acc, WOC = min(f(xg)))
  }))
}

#' Posterior-probability sample size criterion (assurance with separate design and analysis priors)
#'
#' Preposterior probability that \eqn{\Pr(\theta>0\mid X) > \gamma}, normal outcomes, two arms of
#' n. A point-mass design prior (`design_sd = 0`) and flat analysis prior is the ordinary power
#' calculation; a design prior with sd > 0 and flat analysis prior is assurance.
#'
#' @param n Patients per arm (vector).
#' @param design_mean,design_sd Design prior.
#' @param analysis_mean,analysis_sd Analysis prior (`Inf`: flat).
#' @param gamma Posterior probability threshold.
#' @param sigma Outcome sd.
#' @return Data frame with `n` and `prob`.
#' @examples
#' decide_ss_posterior_prob(c(200, 300, 400))
#' @export
decide_ss_posterior_prob <- function(n, design_mean = 0.3, design_sd = 0.1, analysis_mean = 0,
                                     analysis_sd = Inf, gamma = 0.975, sigma = 1) {
  p <- vapply(n, function(ni) {
    a <- .ss_post(ni, sigma, analysis_mean, analysis_sd)
    # success iff w*xbar + (1-w)*ma > qnorm(gamma)*sqrt(pv); xbar ~ N(dm, ds^2 + v)
    thr <- (stats::qnorm(gamma) * sqrt(a$pv) - (1 - a$w) * a$ma) / a$w
    sx <- sqrt(design_sd^2 + a$v)
    1 - stats::pnorm((thr - design_mean) / sx)
  }, numeric(1))
  data.frame(n = n, prob = p)
}

#' Smallest n per arm at which the posterior-probability criterion reaches a target
#'
#' @param target Target probability.
#' @inheritParams decide_ss_posterior_prob
#' @return Smallest integer n per arm.
#' @examples
#' decide_ss_min_n(0.8)
#' @export
decide_ss_min_n <- function(target, ..., nmax = 1e5) {
  # the criterion need not be monotone in n, so scan in chunks for the first n that reaches the target
  nmax <- floor(nmax); start <- 1
  while (start <= nmax) {
    ns <- seq(start, min(start + 9999, nmax))
    hit <- which(decide_ss_posterior_prob(ns, ...)$prob >= target)
    if (length(hit)) return(ns[hit[1]])
    start <- start + 10000
  }
  NA_real_
}

# ---- optimal stopping of a single-arm trial by backward induction ----------------------------

#' Backward induction for a single-arm trial with a beta prior
#'
#' At each n the sponsor stops and adopts when the posterior mean exceeds `p0` (payoff
#' \eqn{M(\hat p - p_0)}) or abandons (0), or pays `c` and enrolls one more patient.
#' Exact operating characteristics are obtained by forward propagation of the count
#' distribution.
#'
#' @param Nmax Maximum patients.
#' @param a,b Beta prior.
#' @param p0 Standard-of-care rate.
#' @param M Value scale of the population.
#' @param c Cost per patient.
#' @param p Optional true response rates for exact operating characteristics.
#' @param n_report Sample sizes at which to report the continuation band.
#' @return List with `table` (data frame `n`, `cont_lo`, `cont_hi`, `adopt_min`), `value`
#'   (value at the start) and `oc` (data frame `p`, `P_adopt`, `E_n`, or NULL).
#' @examples
#' decide_stop_single_arm()$table
#' @export
decide_stop_single_arm <- function(Nmax = 60, a = 3, b = 7, p0 = 0.3, M = 10, c = 0.002,
                                   p = NULL, n_report = c(5, 10, 20, 30, 40, 50, 55)) {
  term <- function(n, x) pmax(M * ((a + x) / (a + b + n) - p0), 0)
  V <- vector("list", Nmax + 1); cont <- vector("list", Nmax + 1)
  V[[Nmax + 1]] <- term(Nmax, 0:Nmax); cont[[Nmax + 1]] <- rep(FALSE, Nmax + 1)
  for (n in (Nmax - 1):0) {
    x <- 0:n; m <- (a + x) / (a + b + n)
    cv <- -c + m * V[[n + 2]][x + 2] + (1 - m) * V[[n + 2]][x + 1]
    t <- term(n, x)
    V[[n + 1]] <- pmax(t, cv); cont[[n + 1]] <- cv > t
  }
  tab <- do.call(rbind, lapply(n_report, function(n) {
    x <- 0:n; cx <- x[cont[[n + 1]]]; ad <- x[term(n, x) > 0]
    data.frame(n = n, cont_lo = if (length(cx)) min(cx) else NA, cont_hi = if (length(cx)) max(cx) else NA,
               adopt_min = if (length(ad)) min(ad) else NA)
  }))
  oc <- NULL
  if (!is.null(p)) {
    oc <- do.call(rbind, lapply(p, function(pp) {
      dist <- 1; Padopt <- 0; En <- 0
      for (n in 0:Nmax) {
        x <- 0:n; stp <- !cont[[n + 1]]
        if (any(stp)) {
          Padopt <- Padopt + sum(dist[stp] * (term(n, x[stp]) > 0))
          En <- En + n * sum(dist[stp])
        }
        if (n == Nmax) break
        d <- dist * !stp  # mass continuing
        dist <- c(d * (1 - pp), 0) + c(0, d * pp)
      }
      data.frame(p = pp, P_adopt = Padopt, E_n = En)
    }))
  }
  list(table = tab, value = V[[1]][1], oc = oc)
}

# ---- pivotal sizing and interim posterior -------------------------------------------------------

#' Pivotal trial sizing: power, assurance and expected NPV
#'
#' Assurance \eqn{A(n)=\Phi((m-z\sqrt{2/n})/\sqrt{2/n+v})} under the posterior \eqn{N(m,v)};
#' \eqn{\mathrm{ENPV}(n) = p_{app} V A(n) - (\text{fixed} + \text{per\_arm}\cdot n)}.
#'
#' @param n Patients per arm.
#' @param m,v Posterior mean and variance of the standardized effect.
#' @param z Critical value.
#' @param fixed Fixed cost.
#' @param per_n Cost per unit of n per arm (two patients at 0.1 each = 0.2, in millions).
#' @param value Value of an approved drug.
#' @param p_approve Probability of approval given success.
#' @return Data frame with `n`, `cost`, `power`, `assurance`, `enpv`.
#' @examples
#' decide_pivotal_table(c(300, 500), 0.3, 0.141^2)
#' @export
decide_pivotal_table <- function(n = c(200, 300, 400, 500, 600, 800, 1000), m, v, z = 1.96,
                                 fixed = 20, per_n = 0.2, value = 1000, p_approve = 0.9) {
  se <- sqrt(2 / n)
  A <- stats::pnorm((m - z * se) / sqrt(se^2 + v))
  data.frame(n = n, cost = fixed + per_n * n, power = stats::pnorm(m / se - z),
             assurance = A, enpv = p_approve * value * A - (fixed + per_n * n))
}

#' ENPV-maximizing pivotal size
#'
#' @inheritParams decide_pivotal_table
#' @param grid Candidate n.
#' @return Named vector `c(n, enpv)` at the grid maximum.
#' @examples
#' decide_pivotal_optimum(0.3, 0.141^2)
#' @export
decide_pivotal_optimum <- function(m, v, grid = 50:2000, ...) {
  t <- decide_pivotal_table(grid, m, v, ...)
  i <- which.max(t$enpv); c(n = t$n[i], enpv = t$enpv[i])
}

#' Interim posterior and predictive probability of final success
#'
#' Normal prior \eqn{N(m_0,s_0^2)} updated by an interim estimate; the predictive probability
#' that the final z statistic (average of two equal-information halves, final se
#' \eqn{se_1/\sqrt2}) exceeds `zfinal`.
#'
#' @param m0,s0 Prior mean and sd (e.g. the skeptical Phase II posterior).
#' @param est,se Interim estimate and se (half information).
#' @param zfinal Final critical value (default 2.02).
#' @return Named list `m`, `s`, `pp`.
#' @examples
#' decide_interim_posterior(0.20, 0.115, 0.15, 0.082)
#' @export
decide_interim_posterior <- function(m0, s0, est, se, zfinal = 2.02) {
  v <- 1 / (1 / s0^2 + 1 / se^2); m <- v * (m0 / s0^2 + est / se^2)
  sef <- se / sqrt(2)
  need <- 2 * zfinal * sef - est  # needed second-half estimate
  list(m = m, s = sqrt(v), pp = stats::pnorm((m - need) / sqrt(v + se^2)))
}

#' Payer's incremental net benefit posterior
#'
#' @param theta,se Posterior mean and sd of the effect.
#' @param qaly_per_unit QALYs per unit effect.
#' @param lambda Willingness to pay per QALY.
#' @param dC Incremental cost.
#' @return List with `mean`, `sd`, `prob_positive`.
#' @examples
#' decide_inb(0.22, 0.058, 0.25, 1e5, 4000)
#' @export
decide_inb <- function(theta, se, qaly_per_unit, lambda, dC) {
  m <- lambda * qaly_per_unit * theta - dC; s <- abs(lambda * qaly_per_unit) * se
  list(mean = m, sd = s, prob_positive = stats::pnorm(m / s))
}

# ---- noninferiority and borrowing worked cases --------------------------------------------------

#' Posterior probability of noninferiority under flat, skeptical and enthusiastic priors
#'
#' Events are bad; \eqn{\theta<M} is noninferior.
#'
#' @param est,se Estimate and standard error.
#' @param M Margin.
#' @param skeptic_sd Sd of the skeptical prior centred at M.
#' @param enthusiast_sd Sd of the enthusiastic prior centred at 0.
#' @return Data frame with `prior`, `post_mean`, `post_sd`, `prob_ni`; the attribute
#'   `upper_limit` holds the frequentist 97.5 percent upper limit.
#' @examples
#' ni_posterior_table(0.01, 0.015, 0.04)
#' @export
ni_posterior_table <- function(est, se, M, skeptic_sd = 0.02, enthusiast_sd = 0.01) {
  upd <- function(pm, ps) {
    if (is.null(ps)) return(c(est, se))
    v <- 1 / (1 / se^2 + 1 / ps^2); c(v * (est / se^2 + pm / ps^2), sqrt(v))
  }
  rows <- list(flat = upd(0, NULL), skeptical = upd(M, skeptic_sd),
               enthusiastic = upd(0, enthusiast_sd))
  out <- data.frame(prior = names(rows), post_mean = vapply(rows, `[`, 1, FUN.VALUE = 0),
                    post_sd = vapply(rows, `[`, 2, FUN.VALUE = 0), row.names = NULL)
  out$prob_ni <- stats::pnorm((M - out$post_mean) / out$post_sd)
  attr(out, "upper_limit") <- est + stats::qnorm(0.975) * se
  out
}

#' Noninferiority margin synthesis
#'
#' With \eqn{\delta_{CP}} the benefit of the control over placebo (positive, from historical
#' trials, \eqn{N(d, s_{CP}^2)}) and \eqn{\theta_{NC}} the new-minus-control event difference
#' (\eqn{N(\text{est}, se^2)}, events bad), the new drug retains at least a fraction f of the
#' control effect iff \eqn{\theta_{NC} < (1-f)\delta_{CP}}. Returns the synthesis probability
#' and the fixed-margin comparison (f times the 95 percent lower limit of \eqn{\delta_{CP}}).
#'
#' @param est,se Current trial estimate and se of \eqn{\theta_{NC}}.
#' @param d,s_cp Historical estimate and se of the control effect over placebo.
#' @param f Fraction to be retained.
#' @return Named list with `prob_retain`, `fixed_margin`, `fixed_ni` (upper 97.5 limit below
#'   fixed margin).
#' @examples
#' ni_synthesis(0.01, 0.015, 0.10, 0.02, 0.5)
#' @export
ni_synthesis <- function(est, se, d, s_cp, f = 0.5) {
  m <- est - (1 - f) * d; s <- sqrt(se^2 + (1 - f)^2 * s_cp^2)
  fm <- (1 - f) * (d - stats::qnorm(0.975) * s_cp)
  list(prob_retain = stats::pnorm(-m / s), fixed_margin = fm,
       fixed_ni = est + stats::qnorm(0.975) * se < fm)
}

#' Rare-disease power-prior table
#'
#' Control posterior \eqn{\mathrm{Beta}(1+a_0y_0+y_C, 1+a_0(n_0-y_0)+n_C-y_C)} (built from
#' [borrow_power_prior()]); treatment \eqn{\mathrm{Beta}(1+y_T,1+n_T-y_T)}. Probabilities by
#' one-dimensional quadrature.
#'
#' @param a0 Discounts.
#' @param yh,nh Historical responders and patients.
#' @param yc,nc Current control; `yt`, `nt` current treatment.
#' @param delta Clinically relevant difference for the second probability.
#' @return Data frame with `a0`, `ess_control`, `p_superior`, `p_delta`.
#' @examples
#' borrow_rare_table()
#' @export
borrow_rare_table <- function(a0 = c(0, 0.25, 0.5, 1), yh = 12, nh = 40, yc = 6, nc = 15,
                              yt = 14, nt = 25, delta = 0.1) {
  pr <- function(ac, bc, d) stats::integrate(function(p)
    stats::dbeta(p, ac, bc) * (1 - stats::pbeta(pmin(p + d, 1), 1 + yt, 1 + nt - yt)), 0, 1,
    rel.tol = 1e-10)$value
  do.call(rbind, lapply(a0, function(a) {
    pp <- borrow_power_prior(yh, nh, a)
    ac <- pp$a + yc; bc <- pp$b + nc - yc
    data.frame(a0 = a, ess_control = a * nh + nc, p_superior = pr(ac, bc, 0),
               p_delta = pr(ac, bc, delta))
  }))
}

#' Log marginal likelihood of current control data as a function of the power-prior discount
#'
#' \eqn{m(a_0) = B(1+a_0y_0+y_C,\,1+a_0(n_0-y_0)+n_C-y_C)/B(1+a_0y_0,\,1+a_0(n_0-y_0))}.
#'
#' @param a0 Grid of discounts.
#' @param yc,nc Current control data.
#' @param yh,nh Historical data.
#' @return Data frame with `a0`, `log_m`, and attribute `a0_hat` (grid maximizer).
#' @examples
#' borrow_marginal_a0(yc = 6)
#' @export
borrow_marginal_a0 <- function(a0 = seq(0, 1, by = 0.01), yc = 6, nc = 15, yh = 12, nh = 40) {
  pp <- borrow_power_prior(yh, nh, a0)
  lm <- lbeta(pp$a + yc, pp$b + nc - yc) - lbeta(pp$a, pp$b)
  out <- data.frame(a0 = a0, log_m = lm)
  attr(out, "a0_hat") <- a0[which.max(lm)]
  out
}

#' Pediatric extrapolation with a mixture prior
#'
#' Prior \eqn{w\,N(\hat\theta_A,s_A^2+\tau^2)+(1-w)N(0,s_v^2)}; observed pediatric estimate
#' \eqn{N(\theta_P,se^2)}.
#'
#' @param thetaA,seA Adult estimate and se.
#' @param tau Allowed adult-pediatric difference sd.
#' @param w Prior weight on extrapolation.
#' @param est,se Pediatric estimate and se.
#' @param sd_vague Sd of the vague component.
#' @return List with `weight` (posterior weight on informative component), `p_benefit`
#'   (\eqn{\Pr(\theta_P<0)}), `p_benefit_alone`.
#' @examples
#' borrow_pediatric_mixture()
#' @export
borrow_pediatric_mixture <- function(thetaA = -0.5, seA = 0.08, tau = 0.15, w = 0.7,
                                     est = -0.35, se = 0.3, sd_vague = 1) {
  s1 <- sqrt(seA^2 + tau^2)
  d1 <- stats::dnorm(est, thetaA, sqrt(s1^2 + se^2)); d0 <- stats::dnorm(est, 0, sqrt(sd_vague^2 + se^2))
  w1 <- w * d1 / (w * d1 + (1 - w) * d0)
  post <- function(m, s) { v <- 1 / (1 / s^2 + 1 / se^2); c(v * (m / s^2 + est / se^2), sqrt(v)) }
  p1 <- post(thetaA, s1); p0 <- post(0, sd_vague)
  list(weight = w1,
       p_benefit = w1 * stats::pnorm(0, p1[1], p1[2]) + (1 - w1) * stats::pnorm(0, p0[1], p0[2]),
       p_benefit_alone = stats::pnorm(-est / se))
}
