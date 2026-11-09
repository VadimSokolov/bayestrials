# Project probabilities onto {q : q >= lower, sum(q) = total} by iterative clipping
# (water-filling): arms below the bound are fixed at it and the remaining arms are
# rescaled proportionally to share what is left, repeated until no arm is below the bound.
rar_bound_probs <- function(p, lower, total = 1) {
  p <- as.matrix(p); J <- ncol(p)
  if (lower * J > total + 1e-12) stop("allocation floor is infeasible for the number of arms")
  base <- pmax(p, 0) + 1e-12
  fixed <- matrix(FALSE, nrow(p), J)
  q <- total * base / rowSums(base)
  for (it in seq_len(J)) {
    low <- !fixed & q < lower - 1e-14
    if (!any(low)) break
    fixed <- fixed | low
    free_mass <- total - lower * rowSums(fixed)
    fb <- base * !fixed
    sc <- fb * (free_mass / pmax(rowSums(fb), 1e-300))
    q <- ifelse(fixed, lower, sc)
  }
  q
}

# Probability that each column of an (reps x M x J) array of posterior draws is the largest.
rar_best_probs <- function(d) {
  reps <- dim(d)[1]; M <- dim(d)[2]; J <- dim(d)[3]
  cur <- matrix(d[, , 1], reps, M); best <- matrix(1L, reps, M)
  if (J > 1L) for (k in 2:J) {
    dk <- matrix(d[, , k], reps, M)
    up <- dk > cur
    best[up] <- k; cur[up] <- dk[up]
  }
  matrix(vapply(seq_len(J), function(k) rowMeans(best == k), numeric(reps)), reps, J)
}

#' Allocation probabilities for the next batch of a response-adaptive trial
#'
#' Batch Thompson sampling with optional burn-in, allocation floor and fixed control share.
#'
#' The "ts" rule allocates each arm with its posterior probability of being the best arm
#' (Beta(1, 1) priors, `M` posterior draws). The "capped" rule uses equal allocation until
#' `burn * N` patients are enrolled and Thompson probabilities afterwards, projected so that every
#' arm receives at least `floor` (two arms) or `floor / K` (more than two arms) while the
#' probabilities still sum to one: arms below the bound are set to it and the others are rescaled
#' proportionally, repeated until no arm is below the bound. For two arms the allocation therefore
#' lies in `[floor, 1 - floor]`. With `fixed_ctrl = TRUE` the control receives exactly `1 / K` and the
#' remaining `1 - 1 / K` is divided among the experimental arms in proportion to their posterior
#' probabilities of being best among the experimental arms (control excluded), with each
#' experimental arm receiving at least `floor / K`.
#'
#' @param rule One of "equal", "ts" (Thompson) or "capped" (burn-in plus floor).
#' @param S,F Matrices (reps x K) of success and failure counts.
#' @param n_done Patients already enrolled.
#' @param N Total sample size.
#' @param K Number of arms (arm 1 is control).
#' @param floor Minimum allocation (divided by K when K > 2) for the capped rule.
#' @param burn Fraction of N allocated equally before adapting (capped rule).
#' @param M Number of posterior draws per replicate.
#' @param fixed_ctrl If TRUE the capped rule fixes the control share at 1/K.
#' @return A reps x K matrix of allocation probabilities.
#' @references Thall, P. F. and Wathen, J. K. (2007). Practical Bayesian adaptive randomisation in clinical trials. European Journal of Cancer 43, 859-866.
#' @examples rar_alloc_probs("ts", matrix(c(3, 6), 1), matrix(c(7, 4), 1), 20, 100, 2)
#' @export
rar_alloc_probs <- function(rule, S, F, n_done, N, K, floor = 0.10, burn = 0.25, M = 200,
                            fixed_ctrl = FALSE) {
  S <- as.matrix(S); F <- as.matrix(F)
  reps <- nrow(S)
  eq <- matrix(1 / K, reps, K)
  if (rule == "equal" || (rule != "ts" && n_done < burn * N)) return(eq)
  lower <- if (K > 2) floor / K else floor
  if (rule == "capped") {
    if (fixed_ctrl && lower * (K - 1) > 1 - 1 / K + 1e-12) stop("allocation floor is infeasible")
    if (!fixed_ctrl && lower * K > 1 + 1e-12) stop("allocation floor is infeasible")
  }
  idx <- rep(seq_len(K), each = M)
  d <- array(stats::rbeta(reps * M * K, S[, idx] + 1, F[, idx] + 1), c(reps, M, K))
  if (rule == "capped" && fixed_ctrl) {
    q <- rar_best_probs(d[, , -1, drop = FALSE])
    q <- rar_bound_probs(q, lower, 1 - 1 / K)
    return(cbind(1 / K, q))
  }
  p <- rar_best_probs(d)
  if (rule == "ts") return(p)
  rar_bound_probs(p, lower, 1)
}

#' Simulate a response-adaptive randomized trial with a binary endpoint
#'
#' Exactly `N` patients are enrolled in batches of `batch`; the last batch is smaller when `N` is
#' not a multiple of `batch`. Time for the drift and the enrolment strata is the batch midpoint as a
#' fraction of `N`.
#'
#' The `rej_best_*` fields refer to the best experimental arm, the experimental arm with the largest
#' true response probability (the lowest-numbered one among ties) compared with control. Under the
#' null (or whenever no experimental arm beats control) these are therefore rejection rates of a
#' true null hypothesis for the first experimental arm. `share_best`, `n_inferior`, `sd_n_best` and
#' `frac_more_on_inferior` refer to the arm with the largest true probability over all arms (control
#' on ties).
#'
#' @param p True response probabilities (arm 1 is control).
#' @param N Sample size.
#' @param rule "equal", "ts" or "capped".
#' @param drift Rise in every arm's response probability from first to last patient.
#' @param reps Number of simulated trials.
#' @param batch Batch size between allocation updates.
#' @param seed Optional seed.
#' @param n_strata Number of enrolment periods for the stratified analysis.
#' @param fixed_ctrl Fix the control share at 1/K (capped rule).
#' @param ... Passed to [rar_alloc_probs()] (floor, burn, M).
#' @return Named list: mean responders (`resp`), share and number on the best arm, sd of n on best,
#'   fraction of trials with fewer than N/K on the best arm, rejection rates of the naive and
#'   period-stratified z tests (best experimental arm, any arm) at one-sided 0.025 Bonferroni, and `share_ctrl`.
#' @references Thall, P., Fox, P. and Wathen, J. (2015). Statistical controversies in clinical research: scientific and ethical problems with adaptive randomization in comparative clinical trials. Annals of Oncology 26, 1621-1628.
#' @examples rar_simulate(c(0.3, 0.45), 100, "ts", reps = 50, seed = 1)$resp
#' @export
rar_simulate <- function(p, N, rule, drift = 0, reps = 4000, batch = 20, seed = NULL,
                         n_strata = 4, fixed_ctrl = FALSE, ...) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(p)
  S <- matrix(0, reps, K); F <- matrix(0, reps, K)
  Ss <- array(0, c(reps, n_strata, K)); Ns <- Ss
  done <- 0L
  while (done < N) {
    size <- as.integer(min(batch, N - done))
    t <- (done + size / 2) / N
    pr <- rar_alloc_probs(rule, S, F, done, N, K, fixed_ctrl = fixed_ctrl, ...)
    cum <- pr
    for (k in seq_len(K)[-1]) cum[, k] <- cum[, k - 1L] + pr[, k]
    u <- matrix(stats::runif(reps * size), reps, size)
    arm <- matrix(1L, reps, size)
    for (k in seq_len(K - 1)) arm <- arm + (u > cum[, k])
    pe <- pmin(pmax(p[arm] + drift * (t - 0.5), 0), 1)
    y <- matrix(as.numeric(stats::runif(reps * size) < pe), reps, size)
    st <- min(as.integer(t * n_strata), n_strata - 1L) + 1L
    for (k in seq_len(K)) {
      m <- arm == k
      S[, k] <- S[, k] + rowSums(y * m); F[, k] <- F[, k] + rowSums((1 - y) * m)
      Ss[, st, k] <- Ss[, st, k] + rowSums(y * m); Ns[, st, k] <- Ns[, st, k] + rowSums(m)
    }
    done <- done + size
  }
  n <- S + F
  zc <- stats::qnorm(1 - 0.025 / (K - 1))
  mx <- function(x) pmax(x, 1)
  zpool <- function(s1, n1, s0, n0) {
    pb <- (s1 + s0) / mx(n1 + n0)
    se <- sqrt(pmax(pb * (1 - pb) * (1 / mx(n1) + 1 / mx(n0)), 1e-12))
    (s1 / mx(n1) - s0 / mx(n0)) / se
  }
  z_naive <- sapply(2:K, function(k) zpool(S[, k], n[, k], S[, 1], n[, 1]))
  z_naive <- matrix(z_naive, reps, K - 1)
  z_strat <- sapply(2:K, function(k) {
    s1 <- matrix(Ss[, , k], reps, n_strata); n1 <- matrix(Ns[, , k], reps, n_strata)
    s0 <- matrix(Ss[, , 1], reps, n_strata); n0 <- matrix(Ns[, , 1], reps, n_strata)
    nk <- n1 + n0
    w <- ifelse(nk > 0, n1 * n0 / mx(nk), 0)
    d <- s1 / mx(n1) - s0 / mx(n0)
    pb <- (s1 + s0) / mx(nk)
    v <- pb * (1 - pb) * (1 / mx(n1) + 1 / mx(n0))
    sw <- pmax(rowSums(w), 1e-12)
    est <- rowSums(w * d) / sw
    vr <- rowSums(w^2 * v) / sw^2
    est / sqrt(pmax(vr, 1e-12))
  })
  z_strat <- matrix(z_strat, reps, K - 1)
  best <- which.max(p)
  best_exp <- which.max(p[-1])  # column of z_*: best experimental arm vs control
  psd <- function(x) sqrt(mean((x - mean(x))^2))
  list(
    resp = mean(rowSums(S)),
    share_best = mean(n[, best] / N),
    n_inferior = mean(N - n[, best]),
    sd_n_best = psd(n[, best]),
    frac_more_on_inferior = mean(n[, best] < N / K),
    rej_best_naive = mean(z_naive[, best_exp] > zc),
    rej_best_strat = mean(z_strat[, best_exp] > zc),
    rej_any_naive = mean(rowSums(z_naive > zc) > 0),
    rej_any_strat = mean(rowSums(z_strat > zc) > 0),
    share_ctrl = mean(n[, 1] / N)
  )
}

rar_rules <- function(K) {
  r <- data.frame(rule_name = c("Equal 1:1", "Thompson", "Capped Thompson"),
                  rule = c("equal", "ts", "capped"), fixed_ctrl = FALSE,
                  stringsAsFactors = FALSE)
  if (K > 2) {
    r$rule_name[1] <- "Equal"
    r <- rbind(r, data.frame(rule_name = "Capped, control fixed at 1/4", rule = "capped",
                             fixed_ctrl = TRUE, stringsAsFactors = FALSE))
  }
  r
}

rar_table <- function(p, N, reps, drift, seed, rules) {
  do.call(rbind, lapply(seq_len(nrow(rules)), function(i) {
    r <- rar_simulate(p, N, rules$rule[i], drift = drift, reps = reps, seed = seed,
                      fixed_ctrl = rules$fixed_ctrl[i])
    data.frame(rule_name = rules$rule_name[i], drift = drift, as.data.frame(r),
               stringsAsFactors = FALSE)
  }))
}

#' Two-arm comparison of equal, Thompson and capped Thompson allocation
#'
#' @param reps Simulated trials.
#' @param N Sample size.
#' @param p True response probabilities (control, experimental).
#' @param seed Optional seed.
#' @return A data.frame with one row per rule.
#' @examples rar_two_arm(reps = 20, N = 100, seed = 1)$resp
#' @export
rar_two_arm <- function(reps = 5000, N = 400, p = c(0.30, 0.45), seed = 11) {
  rar_table(p, N, reps, 0, seed, rar_rules(2))
}

#' Time-drift sweep for response-adaptive allocation under the null
#'
#' @param drifts Drift values.
#' @param reps Simulated trials.
#' @param N Sample size.
#' @param seed Optional seed (reset for each rule and drift).
#' @return A data.frame with one row per rule and drift.
#' @examples rar_drift_sweep(drifts = 0.1, reps = 20, N = 100, seed = 1)$rej_best_naive
#' @export
rar_drift_sweep <- function(drifts = c(0, 0.1, 0.2, 0.3), reps = 5000, N = 400, seed = 12) {
  rules <- rar_rules(2)
  do.call(rbind, lapply(drifts, function(d) rar_table(c(0.30, 0.30), N, reps, d, seed, rules)))
}

#' Four-arm comparison under the alternative and under the null with drift
#'
#' @param reps Simulated trials.
#' @param N Sample size.
#' @param p True response probabilities under the alternative.
#' @param drift Drift applied in the null scenario.
#' @param seed Optional base seed (alternative uses seed, null uses seed + 1).
#' @return A list with data.frames `alt` and `null`, one row per rule.
#' @examples rar_four_arm(reps = 20, N = 100, seed = 1)$alt$resp
#' @export
rar_four_arm <- function(reps = 3000, N = 600, p = c(0.30, 0.30, 0.40, 0.50), drift = 0.2,
                         seed = 13) {
  rules <- rar_rules(4)
  list(alt = rar_table(p, N, reps, 0, seed, rules),
       null = rar_table(rep(0.30, 4), N, reps, drift, if (is.null(seed)) NULL else seed + 1, rules))
}

#' Gittins index for a Bernoulli arm with a beta posterior
#'
#' Bisection on the retirement reward, with dynamic programming over a finite
#' horizon of `D` further patients and geometric discounting.
#'
#' @param a,b Beta posterior parameters (successes + 1, failures + 1 under a uniform prior).
#' @param gamma Discount factor per patient.
#' @param D Truncation horizon of the dynamic program.
#' @param iter Number of bisection steps.
#' @return The Gittins index (a number in (0, 1)).
#' @references Gittins (1979) Bandit processes and dynamic allocation indices. JRSS B 41, 148-177.
#' @examples rar_gittins_index(1, 1, D = 40, iter = 20)
#' @export
rar_gittins_index <- function(a, b, gamma = 0.9, D = 120, iter = 50) {
  lo <- 0; hi <- 1
  for (i in seq_len(iter)) {
    m <- (lo + hi) / 2; Rv <- m / (1 - gamma); V <- rep(Rv, D + 1)
    for (n in (D - 1):0) {
      s <- 0:n; p <- (a + s) / (a + b + n)
      cont <- p * (1 + gamma * V[s + 2]) + (1 - p) * gamma * V[s + 1]
      if (n == 0) c0 <- cont[1]
      V <- pmax(Rv, cont)
    }
    if (c0 > Rv) lo <- m else hi <- m
  }
  (lo + hi) / 2
}

#' Bandit regret simulation: patients on arm 1 over time
#'
#' Simulates `reps` two-arm Bernoulli bandits in parallel and returns the
#' average number of patients assigned to arm 1 after each patient.
#'
#' @param policies Any of "Equal", "UCB1", "Thompson".
#' @param T Horizon (number of patients).
#' @param reps Number of simulated trials.
#' @param p Response probabilities; arm 1 is the inferior arm in the book usage.
#' @param seed Optional seed, set once before all policies are run.
#' @return A T by length(policies) matrix with named columns.
#' @references Thompson (1933); Auer, Cesa-Bianchi and Fischer (2002) Machine Learning 47, 235-256.
#' @examples rar_regret_sim("Thompson", T = 50, reps = 20, seed = 1)[50, ]
#' @export
rar_regret_sim <- function(policies = c("Equal", "UCB1", "Thompson"), T = 5000, reps = 300,
                           p = c(0.30, 0.45), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  one <- function(pol) {
    S <- N <- matrix(0, reps, 2); inf <- numeric(T)
    for (t in seq_len(T)) {
      a <- switch(pol,
        Thompson = max.col(matrix(stats::rbeta(2 * reps, S + 1, N - S + 1), reps)),
        UCB1 = { u <- S / pmax(N, 1) + sqrt(2 * log(t) / pmax(N, 1)); u[N == 0] <- Inf; max.col(u, "first") },
        Equal = sample.int(2, reps, TRUE),
        stop("unknown policy: ", pol))
      i <- cbind(seq_len(reps), a)
      S[i] <- S[i] + (stats::runif(reps) < p[a]); N[i] <- N[i] + 1
      inf[t] <- mean(N[, 1])
    }
    inf
  }
  out <- vapply(policies, one, numeric(T))
  if (T == 1) out <- matrix(out, 1, dimnames = list(NULL, policies))
  out
}

#' Randomized play-the-winner urn (ECMO scenario) versus fixed allocation
#'
#' Runs the urn with one ball per arm, adds a ball for the winning arm after
#' a success (or for the other arm after a failure), then compares the power of a
#' one-sided Fisher exact test under the urn and under fixed equal allocation.
#'
#' @param reps Number of simulated urn trials.
#' @param n Patients per trial.
#' @param p Survival probabilities, c(control, ECMO).
#' @param m Number of trials used for the Fisher power estimates.
#' @param alpha One-sided level.
#' @param seed Optional seed.
#' @return A list with `n_cmt` (mean patients on control), `p_le1` (probability at most one on control), `power_rpw` and `power_fixed`.
#' @references Wei and Durham (1978) JASA 73, 840-843; Bartlett et al. (1985) Pediatrics 76, 479-487.
#' @examples rar_ecmo_urn(reps = 200, m = 50, seed = 1)$n_cmt
#' @export
rar_ecmo_urn <- function(reps = 20000, n = 12, p = c(0.2, 0.8), m = 2000, alpha = 0.025, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  cnt <- s <- matrix(0, reps, 2); balls <- matrix(1, reps, 2)
  for (i in seq_len(n)) {
    a <- 1 + (stats::runif(reps) < balls[, 2] / rowSums(balls)); y <- stats::runif(reps) < p[a]
    j <- cbind(seq_len(reps), a); cnt[j] <- cnt[j] + 1; s[j] <- s[j] + y
    j <- cbind(seq_len(reps), ifelse(y, a, 3 - a)); balls[j] <- balls[j] + 1
  }
  fish <- function(s1, n1, s0, n0)
    stats::fisher.test(matrix(c(s1, n1 - s1, s0, n0 - s0), 2, byrow = TRUE),
                       alternative = "greater")$p.value < alpha
  m <- min(m, reps); ok <- cnt[1:m, 1] > 0 & cnt[1:m, 2] > 0
  rej_rpw <- mean(vapply(1:m, function(i) ok[i] && fish(s[i, 2], cnt[i, 2], s[i, 1], cnt[i, 1]), logical(1)))
  h <- n %/% 2
  ss <- cbind(stats::rbinom(m, h, p[1]), stats::rbinom(m, h, p[2]))
  rej_eq <- mean(vapply(1:m, function(i) fish(ss[i, 2], h, ss[i, 1], h), logical(1)))
  list(n_cmt = mean(cnt[, 1]), p_le1 = mean(cnt[, 1] <= 1), power_rpw = rej_rpw, power_fixed = rej_eq)
}
