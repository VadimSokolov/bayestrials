# E-processes for anytime-valid inference and their simulations.

#' Half-normal mixture e-process
#'
#' Closed-form mixture e-process for H0: theta <= 0 when `S` is the running sum
#' of N(theta, 1) observations. Computed on the log scale
#' ([ev_log_mixture()]) and exponentiated, so it is finite (or correctly `Inf`
#' or 0) wherever the e-value is representable, for example
#' `ev_mixture(-1000, 100, 1)` is about 0.000798.
#'
#' @param S Running sum(s).
#' @param n Sample size(s), nonnegative.
#' @param rho Mixing standard deviation, positive.
#' @param sign `1` for the one-sided (theta > 0) alternative, `-1` for the other side.
#' @return The e-value(s).
#' @examples
#' ev_mixture(5, 100, 0.15)
#' @export
ev_mixture <- function(S, n, rho, sign = 1) {
  exp(ev_log_mixture(S, n, rho, sign))
}

#' Log of the half-normal mixture e-process
#'
#' @inheritParams ev_mixture
#' @return The log e-value(s), numerically stable for large `|S|`.
#' @examples
#' ev_log_mixture(5, 100, 0.15)
#' @export
ev_log_mixture <- function(S, n, rho, sign = 1) {
  if (anyNA(rho) || any(rho <= 0) || any(!is.finite(rho))) stop("`rho` must be positive and finite")
  if (anyNA(n) || any(n < 0)) stop("`n` must be nonnegative")
  if (!all(sign %in% c(-1, 1))) stop("`sign` must be 1 or -1")
  A <- n + 1 / rho^2
  z <- sign * S / sqrt(A)
  # log Phi(z) + z^2 / 2 without forming Inf - Inf for large |z|: for z << 0 use the Mills-ratio form
  zz <- pmin(z, -30)
  lm <- ifelse(z < -30,
               -log(-zz) - 0.5 * log(2 * pi) + log1p(-1 / zz^2 + 3 / zz^4 - 15 / zz^6 + 105 / zz^8),
               stats::pnorm(z, log.p = TRUE) + z^2 / 2)
  log(2 / (rho * sqrt(A))) + lm
}

#' Random walks of N(delta, 1) increments
#'
#' @param delta Per-observation mean.
#' @param nmax Walk length.
#' @param reps Number of walks.
#' @param seed Optional seed.
#' @return A `reps` by `nmax` matrix of cumulative sums (one walk per row).
#' @examples
#' dim(ev_walks(0, 10, 5))
#' @export
ev_walks <- function(delta, nmax, reps, seed = NULL) {
  stopifnot(length(nmax) == 1L, nmax >= 1, nmax == round(nmax),
            length(reps) == 1L, reps >= 1, reps == round(reps))
  if (!is.null(seed)) set.seed(seed)
  M <- matrix(stats::rnorm(reps * nmax, delta, 1), reps, nmax)
  for (j in seq_len(nmax)[-1]) M[, j] <- M[, j] + M[, j - 1]
  M
}

#' First crossing index
#'
#' @param mask Logical matrix, one row per path.
#' @return For each row, the first column index where `mask` is `TRUE`, or
#'   `ncol(mask) + 1` if never.
#' @examples
#' ev_first_cross(matrix(c(FALSE, TRUE, FALSE, FALSE), 2, byrow = TRUE))
#' @export
ev_first_cross <- function(mask) {
  nmax <- ncol(mask)
  ifelse(rowSums(mask) > 0, max.col(mask * 1, ties.method = "first"), nmax + 1)
}

ev_chunks <- function(reps, chunk) {
  starts <- seq(1, reps, by = chunk)
  pmin(chunk, reps - starts + 1)
}

#' Type I error of repeated z tests versus the e-process
#'
#' Under the null, the false positive rate of one-sided z tests at level
#' `alpha` at `k` equally spaced looks, and of the mixture e-process
#' crossing `1/alpha` at any time.
#'
#' @param nmax Maximum sample size.
#' @param looks_list Vector of numbers of looks.
#' @param reps Number of simulated paths.
#' @param rho Mixing standard deviation.
#' @param alpha One-sided level.
#' @param seed Optional seed.
#' @return A list: `naive` (data.frame with `looks`, `false_positive`) and
#'   `evalue` (e-process false positive rate).
#' @examples
#' ev_naive_false_alarm(50, c(1, 5), reps = 200, seed = 1)$naive
#' @export
ev_naive_false_alarm <- function(nmax, looks_list, reps = 20000, rho = 0.15,
                                 alpha = 0.025, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(nmax) == 1L, nmax >= 1, nmax == round(nmax), alpha > 0, alpha < 1,
            length(rho) == 1L, rho > 0, reps >= 1,
            all(looks_list >= 1), all(looks_list == round(looks_list)), all(looks_list <= nmax))
  n <- seq_len(nmax)
  zc <- stats::qnorm(1 - alpha)
  hits <- numeric(length(looks_list)); e_hit <- 0
  for (m in ev_chunks(reps, 4000)) {
    S <- ev_walks(0, nmax, m)
    Z <- sweep(S, 2, sqrt(n), "/")
    for (i in seq_along(looks_list)) {
      k <- looks_list[i]
      idx <- unique(round(seq(nmax / k, nmax, length.out = k)))
      hits[i] <- hits[i] + sum(rowSums(Z[, idx, drop = FALSE] >= zc) > 0)
    }
    LE <- ev_log_mixture(S, matrix(n, m, nmax, byrow = TRUE), rho)
    e_hit <- e_hit + sum(rowSums(LE >= log(1 / alpha)) > 0)
  }
  list(naive = data.frame(looks = looks_list, false_positive = hits / reps),
       evalue = e_hit / reps)
}

#' SELECT-like comparison of Bayes, e-value and O'Brien-Fleming rules
#'
#' Event-time Gaussian approximation: per-event drift `-log(hr)/2`, prior sd
#' `tau/2` on that scale. The Bayes posterior cutoff and the four-look OBF
#' constant are calibrated by simulation on the null.
#'
#' @param nmax Maximum number of events.
#' @param hr Hazard ratio under the alternative.
#' @param tau Prior sd on the log hazard ratio scale.
#' @param reps Number of paths per scenario.
#' @param alpha One-sided level.
#' @param seed Optional seed.
#' @return A list: `cb` (z cutoff for the posterior of delta, so the posterior
#'   probability cutoff is `pnorm(cb)`), `co` (OBF constant), `nmax`, and
#'   `null` and `alt`, each a list of vectors of stopping events (`nmax + 1`
#'   if no action) for `bayes`, `evalue`, `obf`.
#' @examples
#' ev_select_like(nmax = 100, reps = 200, seed = 1)$co
#' @export
ev_select_like <- function(nmax = 1270, hr = 0.80, tau = 0.25, reps = 20000,
                           alpha = 0.025, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(length(nmax) == 1L, nmax >= 4, nmax == round(nmax), hr > 0, tau > 0,
            alpha > 0, alpha < 1, reps >= 1)
  n <- seq_len(nmax)
  delta <- -log(hr) / 2
  rho <- tau / 2
  looks <- c(0.25, 0.5, 0.75, 1)
  li <- as.integer(looks * nmax)
  zpost <- function(S) sweep(S, 2, sqrt(1 + n * rho^2), "/") * rho
  zlook <- function(S) sweep(S[, li, drop = FALSE], 2, sqrt(n[li]) / sqrt(looks), "/")
  chunks <- ev_chunks(reps, 2000)
  mb <- mo <- NULL
  for (m in chunks) {
    S <- ev_walks(0, nmax, m)
    mb <- c(mb, do.call(pmax, as.data.frame(zpost(S))))
    mo <- c(mo, do.call(pmax, as.data.frame(zlook(S))))
  }
  cb <- stats::quantile(mb, 1 - alpha, names = FALSE)
  co <- stats::quantile(mo, 1 - alpha, names = FALSE)
  lE <- log(1 / alpha)
  res <- list()
  for (nm in c("null", "alt")) {
    d <- if (nm == "null") 0 else delta
    st <- list(bayes = NULL, evalue = NULL, obf = NULL)
    for (m in chunks) {
      S <- ev_walks(d, nmax, m)
      nb <- ev_first_cross(zpost(S) >= cb)
      LE <- ev_log_mixture(S, matrix(n, m, nmax, byrow = TRUE), rho)
      ne <- ev_first_cross(LE >= lE)
      ko <- ev_first_cross(zlook(S) >= co)
      no <- ifelse(ko <= 4, n[li][pmin(ko, 4)], nmax + 1)
      st$bayes <- c(st$bayes, nb); st$evalue <- c(st$evalue, ne); st$obf <- c(st$obf, no)
    }
    res[[nm]] <- st
  }
  list(cb = cb, co = co, nmax = nmax, null = res$null, alt = res$alt)
}

#' Betting martingale log wealth paths
#'
#' Log wealth of the martingale \eqn{\prod(1 + \lambda D_i)} with
#' \eqn{D = } outcome(T) minus outcome(C) for Bernoulli pairs.
#'
#' @param pT,pC Response rates in treatment and control.
#' @param lam Betting fraction in `[-1, 1]`: each increment is -1, 0 or 1, so
#'   `1 + lam * D >= 0` (nonnegative wealth) exactly when `abs(lam) <= 1`.
#' @param npairs Number of pairs.
#' @param reps Number of paths.
#' @param seed Optional seed.
#' @return A `reps` by `npairs` matrix of log wealth.
#' @examples
#' dim(ev_betting_paths(0.45, 0.30, 0.31, 20, 5, seed = 1))
#' @export
ev_betting_paths <- function(pT, pC, lam, npairs, reps, seed = NULL) {
  if (!all(c(pT, pC) >= 0 & c(pT, pC) <= 1)) stop("`pT` and `pC` must be probabilities in [0, 1]")
  if (length(lam) != 1L || !is.finite(lam) || abs(lam) > 1)
    stop("`lam` must be a single number in [-1, 1] (nonnegative wealth)")
  stopifnot(length(npairs) == 1L, npairs >= 1, npairs == round(npairs),
            length(reps) == 1L, reps >= 1, reps == round(reps))
  if (!is.null(seed)) set.seed(seed)
  D <- (matrix(stats::runif(reps * npairs), reps, npairs) < pT) * 1 -
       (matrix(stats::runif(reps * npairs), reps, npairs) < pC) * 1
  W <- log1p(lam * D)
  for (j in seq_len(npairs)[-1]) W[, j] <- W[, j] + W[, j - 1]
  W
}
