# Simon two-stage design and the MCP-Mod multiple contrast test (helpers used in chapters 13 and 15).

#' Exact operating characteristics of a Simon two-stage design
#'
#' Stage 1 enrolls `n1` patients and stops for futility if at most `r1` respond; otherwise
#' the trial continues to `n` patients and the drug is declared promising if more than `r`
#' respond. Both quantities are exact binomial sums.
#'
#' @param n1,r1,n,r Design: first-stage size and futility bound, total size and final bound.
#' @param p True response rate.
#' @return Named vector `c(P_reject, EN)`: probability of declaring the drug promising and
#'   expected sample size at response rate `p`.
#' @export
#' @examples
#' simon_oc(13, 3, 43, 12, 0.2)
simon_oc <- function(n1, r1, n, r, p) {
  x1 <- (r1 + 1):n1
  c(P_reject = sum(stats::dbinom(x1, n1, p) * stats::pbinom(r - x1, n - n1, p, lower.tail = FALSE)),
    EN = n1 + (1 - stats::pbinom(r1, n1, p)) * (n - n1))
}

#' Simon's optimal two-stage design
#'
#' Exact search over `(n1, r1, n, r)` for the design with type I error at most `alpha` at
#' `p0` and power at least `power` at `p1` that has the smallest expected sample size under `p0`.
#'
#' @param p0,p1 Uninteresting and target response rates.
#' @param alpha Maximum type I error (at `p0`).
#' @param power Minimum power (at `p1`).
#' @param nmax Largest total sample size searched.
#' @return One-row data frame with `n1`, `r1`, `n`, `r` and `EN`, the expected sample size at `p0`.
#' @export
#' @examples
#' simon_optimal(0.2, 0.4, 0.05, 0.80)
simon_optimal <- function(p0, p1, alpha, power, nmax = 50) {
  best <- NULL
  for (n1 in 5:(nmax - 1)) for (r1 in 0:min(n1 - 1, floor(n1 * p1))) for (n in (n1 + 1):nmax) {
    r <- r1:min(n, n1 + r1 + 20)
    a <- vapply(r, function(k) simon_oc(n1, r1, n, k, p0)[[1]], numeric(1))
    i <- which(a <= alpha)[1]
    if (is.na(i) || simon_oc(n1, r1, n, r[i], p1)[[1]] < power) next
    en <- simon_oc(n1, r1, n, r[i], p0)[[2]]
    if (is.null(best) || en < best$EN - 1e-9) best <- data.frame(n1 = n1, r1 = r1, n = n, r = r[i], EN = en)
  }
  best
}

#' Optimal contrasts and critical value of the MCP-Mod multiple contrast test
#'
#' For each candidate dose-response shape (a column of `shapes`, evaluated at the doses) the
#' optimal contrast with equal group sizes is the centered, unit-length shape. The critical
#' value is the one-sided `1 - alpha` quantile of the maximum of the contrast statistics,
#' whose joint null distribution is multivariate normal with correlation `crossprod(contrasts)`
#' (known standard deviation). Needs the `mvtnorm` package.
#'
#' @param shapes Matrix with one column per candidate shape and one row per dose.
#' @param alpha One-sided familywise type I error.
#' @return List with `contrasts` (matrix, same layout as `shapes`) and `crit`.
#' @export
#' @examples
#' d <- c(0, 10, 25, 50, 100)
#' mcp_contrasts(cbind(linear = d / 100, emax = d / (d + 15)))$crit
mcp_contrasts <- function(shapes, alpha = 0.025) {
  if (!requireNamespace("mvtnorm", quietly = TRUE)) stop("package 'mvtnorm' is required")
  ctr <- apply(shapes, 2, function(m) (m - mean(m)) / sqrt(sum((m - mean(m))^2)))
  crit <- mvtnorm::qmvnorm(1 - alpha, tail = "lower.tail", sigma = crossprod(ctr))$quantile
  list(contrasts = ctr, crit = crit)
}

#' Rejection probabilities of the MCP-Mod contrast test and two traditional tests
#'
#' Simulates dose means with known standard deviation and reports the probability of a
#' one-sided signal for the multiple contrast test, for the top dose against the first
#' (placebo) at level `alpha`, and for the best of all pairwise comparisons with placebo
#' with a Bonferroni correction.
#'
#' @param mu True mean at each dose (first entry is placebo).
#' @param contrasts,crit Output of [mcp_contrasts()].
#' @param n Patients per dose.
#' @param B Number of simulated trials.
#' @param sd Known standard deviation.
#' @param alpha One-sided level of the traditional tests.
#' @param seed Optional seed.
#' @return Named vector `c(mcp, top_dose, bonferroni)`.
#' @export
#' @examples
#' d <- c(0, 10, 25, 50, 100)
#' m <- mcp_contrasts(cbind(linear = d / 100, emax = d / (d + 15)))
#' mcp_power(0.5 * d / 100, m$contrasts, m$crit, n = 40, B = 2000, seed = 1)
mcp_power <- function(mu, contrasts, crit, n, B = 20000, sd = 1, alpha = 0.025, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  k <- length(mu)
  ybar <- matrix(stats::rnorm(B * k, rep(mu, each = B), sd / sqrt(n)), B)
  zp <- (ybar[, -1, drop = FALSE] - ybar[, 1]) / (sd * sqrt(2 / n))
  c(mcp = mean(apply(ybar %*% contrasts * sqrt(n) / sd, 1, max) > crit),
    top_dose = mean(zp[, k - 1] > stats::qnorm(1 - alpha)),
    bonferroni = mean(apply(zp, 1, max) > stats::qnorm(1 - alpha / (k - 1))))
}
