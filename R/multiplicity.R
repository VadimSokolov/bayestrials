#' Bonferroni procedure
#'
#' @param p numeric vector of one-sided p-values.
#' @param alpha one-sided familywise level.
#' @return logical vector, `TRUE` where the hypothesis is rejected.
#' @examples mult_bonferroni(c(0.001, 0.02, 0.4))
#' @references Bonferroni (1936).
#' @export
mult_bonferroni <- function(p, alpha = 0.025) {
  p <= alpha / length(p)
}

#' Holm step-down procedure
#'
#' @inheritParams mult_bonferroni
#' @return logical vector of rejections.
#' @examples mult_holm(c(0.001, 0.02, 0.4))
#' @references Holm (1979), Scandinavian Journal of Statistics 6, 65-70.
#' @export
mult_holm <- function(p, alpha = 0.025) {
  m <- length(p)
  o <- order(p)
  rej <- rep(FALSE, m)
  for (r in seq_len(m)) {
    i <- o[r]
    if (p[i] <= alpha / (m - r + 1)) rej[i] <- TRUE else break
  }
  rej
}

#' Hochberg step-up procedure
#'
#' @inheritParams mult_bonferroni
#' @return logical vector of rejections.
#' @examples mult_hochberg(c(0.001, 0.02, 0.4))
#' @references Hochberg (1988), Biometrika 75, 800-802.
#' @export
mult_hochberg <- function(p, alpha = 0.025) {
  m <- length(p)
  o <- order(p, decreasing = TRUE)
  rej <- rep(FALSE, m)
  for (r in seq_len(m)) {
    i <- o[r]
    if (p[i] <= alpha / r) {
      rej[p <= p[i]] <- TRUE
      break
    }
  }
  rej
}

#' Sequentially rejective graphical procedure
#'
#' @param p numeric vector of p-values.
#' @param w initial weights (sum at most 1).
#' @param G transition matrix (rows are the source hypothesis).
#' @param alpha one-sided familywise level.
#' @return logical vector of rejections.
#' @examples mult_graph(c(0.01, 0.02), c(0.5, 0.5), matrix(c(0, 1, 1, 0), 2))
#' @references Bretz, Maurer, Brannath and Posch (2009), Statistics in Medicine 28, 586-604.
#' @export
mult_graph <- function(p, w, G, alpha = 0.025) {
  m <- length(p)
  active <- rep(TRUE, m)
  rej <- rep(FALSE, m)
  G <- as.matrix(G)
  repeat {
    cand <- which(active & p <= w * alpha)
    if (length(cand) == 0) return(rej)
    i <- cand[1]
    rej[i] <- TRUE
    active[i] <- FALSE
    w2 <- w
    G2 <- G
    for (l in seq_len(m)) {
      if (l == i) next
      w2[l] <- w[l] + w[i] * G[i, l]
      for (k in seq_len(m)) {
        if (k == i || k == l) {
          G2[l, k] <- 0
        } else {
          d <- 1 - G[l, i] * G[i, l]
          G2[l, k] <- if (d > 0) (G[l, k] + G[l, i] * G[i, k]) / d else 0
        }
      }
    }
    w2[i] <- 0
    G2[i, ] <- 0
    G2[, i] <- 0
    w <- w2
    G <- G2
  }
}

#' Simulate multiplicity procedures against a shared control
#'
#' `k` one-sided z tests; equicorrelated 0.5 when `rho` is nonzero.
#'
#' @param k number of hypotheses.
#' @param effects effects in sd units (length `k`, recycled).
#' @param n patients per arm.
#' @param alpha one-sided level.
#' @param nsim number of simulated trials.
#' @param rho 0 for independent tests, nonzero for shared-control correlation 0.5.
#' @param seed optional seed.
#' @return data.frame with one row per procedure: `procedure`, `fwer` (P(any false rejection)) and `power` (mean per-hypothesis power among nonnull hypotheses).
#' @examples mult_simulate(3, c(0, 0, 0.3), nsim = 500, seed = 1)
#' @export
mult_simulate <- function(k, effects, n = 100, alpha = 0.025, nsim = 20000, rho = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  eff <- rep_len(effects, k)
  zc <- matrix(rnorm(nsim), nsim, k)
  ze <- matrix(rnorm(nsim * k), nsim, k)
  z <- if (rho != 0) sqrt(0.5) * zc + sqrt(0.5) * ze else ze
  z <- z + matrix(eff * sqrt(n / 2), nsim, k, byrow = TRUE)
  p <- 1 - pnorm(z)
  null <- eff == 0
  procs <- list(unadjusted = function(q) q <= alpha,
                Bonferroni = function(q) mult_bonferroni(q, alpha),
                Holm = function(q) mult_holm(q, alpha),
                Hochberg = function(q) mult_hochberg(q, alpha))
  rows <- lapply(names(procs), function(nm) {
    R <- if (nm == "unadjusted") p <= alpha else t(apply(p, 1, procs[[nm]]))
    R <- matrix(R, nsim, k)
    data.frame(procedure = nm,
               fwer = if (any(null)) mean(rowSums(R[, null, drop = FALSE]) > 0) else NA_real_,
               power = if (any(!null)) mean(R[, !null, drop = FALSE]) else NA_real_)
  })
  do.call(rbind, rows)
}

#' Log odds and variances from the ulcer data
#'
#' @param data data.frame with columns `n` and `recurrences` (default [ulcer]).
#' @return list with `y` (log odds of recurrence) and `v` (variances).
#' @examples shrink_logodds()
#' @export
shrink_logodds <- function(data = ulcer) {
  n <- data$n
  r <- data$recurrences
  list(y = log(r / (n - r)), v = 1 / r + 1 / (n - r))
}

#' Empirical Bayes (marginal likelihood) fit of a normal hierarchy
#'
#' Maximum marginal likelihood estimate of tau for `y_i ~ N(mu, v_i + tau^2)`.
#'
#' @param y estimates.
#' @param v sampling variances.
#' @return list with `mu`, `tau`, `shrunk` (posterior means) and `sd` (posterior sds).
#' @examples l <- shrink_logodds(); shrink_eb_fit(l$y, l$v)$tau
#' @export
shrink_eb_fit <- function(y, v) {
  nll <- function(t) {
    s <- v + t^2
    w <- 1 / s
    mu <- sum(w * y) / sum(w)
    0.5 * sum(log(s) + (y - mu)^2 / s)
  }
  t <- optimize(nll, c(1e-6, 3), tol = 1e-8)$minimum
  w <- 1 / (v + t^2)
  mu <- sum(w * y) / sum(w)
  list(mu = mu, tau = t,
       shrunk = (t^2 * y + v * mu) / (t^2 + v),
       sd = sqrt(1 / (1 / v + 1 / t^2)))
}

#' Posterior probability that a normal effect is positive
#'
#' @param est estimate.
#' @param se standard error.
#' @param prior_sd sd of a mean-zero normal prior.
#' @return list with `prob` (Pr(theta > 0 | data)) and `precision`.
#' @examples shrink_post_prob(1, 0.4, 1)
#' @export
shrink_post_prob <- function(est, se, prior_sd) {
  A <- 1 / se^2 + 1 / prior_sd^2
  list(prob = pnorm((est / se^2) / sqrt(A)), precision = A)
}
