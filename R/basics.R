#' Simulate repeated two-arm binary trials
#'
#' Reruns a whole two-arm trial `nsim` times and returns the observed counts and risk differences.
#'
#' @param p0,p1 Response rates in control and treatment.
#' @param n Patients per arm.
#' @param nsim Number of simulated trials.
#' @param seed Optional seed.
#' @return A data.frame with columns `x0`, `x1` (responders) and `diff` (risk difference).
#' @examples
#' head(sim_two_arm_binary(0.30, 0.45, 100, nsim = 5, seed = 1))
#' @export
sim_two_arm_binary <- function(p0, p1, n, nsim, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x0 <- stats::rbinom(nsim, n, p0)
  x1 <- stats::rbinom(nsim, n, p1)
  data.frame(x0 = x0, x1 = x1, diff = (x1 - x0) / n)
}

#' Pooled two-proportion z statistic
#'
#' Vectorised over `x0` and `x1`. A zero standard error gives a statistic of 0.
#'
#' @param x0,x1 Responders in control and treatment (scalars or vectors).
#' @param n Patients per arm.
#' @return Numeric vector of z statistics.
#' @examples
#' z_test_binary(30, 45, 100)
#' @export
z_test_binary <- function(x0, x1, n) {
  pbar <- (x0 + x1) / (2 * n)
  se <- sqrt(2 * pbar * (1 - pbar) / n)
  se[se == 0] <- Inf
  ((x1 - x0) / n) / se
}

#' Wald confidence interval for a risk difference
#'
#' @param x0,x1 Responders in control and treatment (scalars or vectors).
#' @param n Patients per arm.
#' @param level Confidence level.
#' @return A data.frame with columns `lower` and `upper`.
#' @examples
#' wald_ci(30, 45, 100)
#' @export
wald_ci <- function(x0, x1, n, level = 0.95) {
  d <- (x1 - x0) / n
  p0 <- x0 / n
  p1 <- x1 / n
  se <- sqrt(p0 * (1 - p0) / n + p1 * (1 - p1) / n)
  z <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(lower = d - z * se, upper = d + z * se)
}

#' Randomization test for a difference in means
#'
#' Shuffles the arm labels `nsim` times and recomputes the difference in means.
#'
#' @param y1,y0 Outcomes in treatment and control.
#' @param nsim Number of shuffles.
#' @param seed Optional seed.
#' @return A list with `obs`, `sims` (shuffled differences), `p_one` and `p_two`.
#' @examples
#' perm_test(c(5, 6, 7, 8), c(1, 2, 3, 4), nsim = 200, seed = 1)$p_one
#' @export
perm_test <- function(y1, y0, nsim = 20000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  y <- c(y1, y0)
  n1 <- length(y1)
  n0 <- length(y0)
  tot <- sum(y)
  obs <- mean(y1) - mean(y0)
  s1 <- vapply(seq_len(nsim), function(i) sum(y[sample.int(length(y), n1)]), numeric(1))
  sims <- s1 / n1 - (tot - s1) / n0
  list(obs = obs, sims = sims,
       p_one = (sum(sims >= obs - 1e-12) + 1) / (nsim + 1),
       p_two = (sum(abs(sims) >= abs(obs) - 1e-12) + 1) / (nsim + 1))
}
