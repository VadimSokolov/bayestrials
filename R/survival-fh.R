#' Fleming-Harrington weighted log-rank statistics for two groups
#'
#' Computes the standardized weighted log-rank statistics with Fleming-Harrington
#' weights \eqn{S(t)^\rho \{1 - S(t)\}^\gamma}, where \eqn{S} is the pooled
#' Kaplan-Meier estimate just before each event time, for the four pairs
#' \eqn{(\rho, \gamma) \in \{(0,0), (0,1), (1,0), (1,1)\}}, together with the
#' correlation matrix of the four statistics. The first statistic is the ordinary
#' log-rank test. The maximum of the four, referred to the multivariate normal
#' distribution with this correlation matrix, is the MaxCombo test.
#'
#' Each statistic is the weighted sum of observed minus expected events in group 0,
#' divided by the square root of its hypergeometric variance. A positive value
#' means more events than expected in group 0, that is, an advantage for group 1.
#'
#' @param time Follow-up times.
#' @param status Event indicator (1 = event, 0 = censored).
#' @param grp Group indicator (0 or 1).
#' @return A list with `z`, the four standardized statistics in the order
#'   (0,0), (0,1), (1,0), (1,1), and `corr`, their 4 by 4 correlation matrix.
#' @references Harrington DP, Fleming TR (1982). A class of rank test procedures
#'   for censored survival data. Biometrika 69:553-566.
#'   Lee S-H (2007). On the versatility of the combination of the weighted
#'   log-rank statistics. Computational Statistics and Data Analysis 51:6557-6564.
#' @examples
#' set.seed(1)
#' tm <- c(rexp(100, 0.15), rexp(100, 0.1))
#' fh_z(tm, rep(1, 200), rep(0:1, each = 100))$z
#' @export
fh_z <- function(time, status, grp) {
  ut <- sort(unique(time[status == 1]))
  d  <- tabulate(match(time[status == 1], ut), length(ut))
  d0 <- tabulate(match(time[status == 1 & grp == 0], ut), length(ut))
  n  <- length(time) - findInterval(ut, sort(time), left.open = TRUE)
  n0 <- sum(grp == 0) - findInterval(ut, sort(time[grp == 0]), left.open = TRUE)
  S  <- c(1, cumprod(1 - d / n))[seq_along(ut)]
  v  <- d * (n0 / n) * (1 - n0 / n) * ifelse(n > 1, (n - d) / (n - 1), 0)
  W  <- cbind(1, 1 - S, S, S * (1 - S))        # FH(0,0), (0,1), (1,0), (1,1)
  V  <- crossprod(W, W * v)
  list(z = colSums(W * (d0 - d * n0 / n)) / sqrt(diag(V)), corr = stats::cov2cor(V))
}
