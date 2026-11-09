# Predictive probability, conditional power, Lee-Liu monitoring and backward induction.
# Port of code/ch13_predictive.py.

.pp_cache <- new.env(parent = emptyenv())

# Pr(p1 > p0 + m) on the full grid of states, by the rank-scale integral of
# prob_superior(method = "exact"): rows index treatment Beta(A1, B1), columns control
# Beta(A0, B0). When both arms share the same states and shapes and m = 0 (exchangeable),
# only the lower triangle is integrated and the rest follows from Pr(p1 > p0) = 1 - Pr(p0 > p1)
# (the diagonal is exactly 1/2). Rows are computed in parallel (forked workers, Unix only) when
# `options(bayestrials.cores = k)` with k > 1; the default is serial. Results are cached for the
# session (`cache = FALSE` skips it).
.pr_superior_grid <- function(A1, B1, A0, B0, m = 0, cache = TRUE) {
  key <- paste(sprintf("%.17g", c(length(A1), A1, B1, A0, B0, m)), collapse = "|")
  if (cache && !is.null(hit <- .pp_cache[[key]])) return(hit)
  n1 <- length(A1); n0 <- length(A0)
  out <- matrix(0, n1, n0)
  sym <- m == 0 && n1 == n0 && identical(as.numeric(A1), as.numeric(A0)) &&
    identical(as.numeric(B1), as.numeric(B0))
  row <- function(i) {
    r <- numeric(n0)
    jj <- if (sym) seq_len(i - 1L) else seq_len(n0)
    for (j in jj) r[j] <- prob_superior_exact(A1[i], B1[i], A0[j], B0[j], m)
    r
  }
  cores <- getOption("bayestrials.cores", 1L)
  rows <- if (cores > 1L && .Platform$OS.type == "unix" && n1 > 20L) {
    parallel::mclapply(seq_len(n1), row, mc.cores = cores)
  } else lapply(seq_len(n1), row)
  for (i in seq_len(n1)) out[i, ] <- rows[[i]]
  if (sym) {
    out[upper.tri(out)] <- 1 - t(out)[upper.tri(out)]
    diag(out) <- 0.5
  }
  if (cache) assign(key, out, envir = .pp_cache)
  out
}

#' Conditional power on the normal (B-value) scale
#'
#' `P(Z_1 > za | Z_t = z)` when the drift (`E Z_1`) is `theta`.
#'
#' @param z Current z statistic.
#' @param t Information fraction in (0, 1).
#' @param theta Assumed drift.
#' @param za Critical value (default `qnorm(0.975)`).
#' @return Numeric vector of probabilities.
#' @export
#' @references Lan, K. K. G., Simon, R., Halperin, M. (1982). Stochastically curtailed tests in
#'   long-term clinical trials. Sequential Analysis 1, 207-219.
#' @examples
#' cond_power(1.5, 0.5, 2)
cond_power <- function(z, t, theta, za = stats::qnorm(0.975)) {
  b <- z * sqrt(t)
  stats::pnorm((za - b - theta * (1 - t)) / sqrt(1 - t), lower.tail = FALSE)
}

#' Predictive probability of success on the normal scale
#'
#' Averages [cond_power()] over the posterior of the drift `theta`, prior `N(m, s2)`
#' (`s2 = Inf` is flat), after observing `B = z sqrt(t)`.
#'
#' @param z Current z statistic.
#' @param t Information fraction.
#' @param m,s2 Prior mean and variance of the drift.
#' @param za Critical value.
#' @return Numeric vector of probabilities.
#' @export
#' @references Spiegelhalter, D. J., Freedman, L. S., Blackburn, P. R. (1986). Monitoring
#'   clinical trials: conditional or predictive power? Controlled Clinical Trials 7, 8-17.
#' @examples
#' pp_normal(1.5, 0.5)
pp_normal <- function(z, t, m = 0, s2 = Inf, za = stats::qnorm(0.975)) {
  b <- z * sqrt(t)
  prec0 <- if (is.infinite(s2)) 0 else 1 / s2
  prior_term <- if (is.infinite(s2)) 0 else m / s2
  vp <- 1 / (t + prec0)
  mp <- (b + prior_term) * vp
  mean <- b + (1 - t) * mp
  v <- (1 - t) + (1 - t)^2 * vp
  stats::pnorm((za - mean) / sqrt(v), lower.tail = FALSE)
}

#' Exact predictive probability, one arm
#'
#' Probability that, at `N` patients, `Pr(p > p0 | data) > eta`, given `s` responses in `n`.
#'
#' @param s,n Responses and patients so far.
#' @param N Final sample size.
#' @param p0 Null response rate.
#' @param eta Posterior probability threshold.
#' @param a,b Beta prior.
#' @return A probability.
#' @export
#' @examples
#' pp_one_arm(8, 20, 40, 0.2, 0.95)
pp_one_arm <- function(s, n, N, p0, eta, a = 1, b = 1) {
  m <- N - n
  x <- 0:m
  w <- .dbetabinom(x, m, a + s, b + n - s)
  post <- stats::pbeta(p0, a + s + x, b + N - s - x, lower.tail = FALSE)
  sum(w * (post > eta))
}

#' Lee-Liu monitoring table
#'
#' For each look `n`, the largest `s` that stops for futility (predictive probability below
#' `th_fut`) and the smallest `s` that stops for efficacy (above `th_eff`).
#'
#' @param N Final sample size.
#' @param p0 Null response rate.
#' @param looks Integer vector of look sample sizes.
#' @param th_fut,th_eff Futility and efficacy thresholds on the predictive probability.
#' @param eta Final posterior probability threshold.
#' @param a,b Beta prior.
#' @return Data frame with columns `n`, `s_fut`, `s_eff` (`NA` when no such `s`).
#' @export
#' @references Lee, J. J., Liu, D. D. (2008). A predictive probability design for phase II
#'   cancer clinical trials. Clinical Trials 5, 93-106.
#' @examples
#' lee_liu_table(20, 0.2, c(10, 15), 0.05, 0.95, 0.9)
lee_liu_table <- function(N, p0, looks, th_fut, th_eff, eta, a = 1, b = 1) {
  res <- lapply(looks, function(n) {
    pp <- vapply(0:n, function(s) pp_one_arm(s, n, N, p0, eta, a, b), numeric(1))
    s <- 0:n
    fut <- s[pp < th_fut]; eff <- s[pp > th_eff]
    data.frame(n = n, s_fut = if (length(fut)) max(fut) else NA_integer_,
               s_eff = if (length(eff)) min(eff) else NA_integer_)
  })
  do.call(rbind, res)
}

#' Exact operating characteristics of a Lee-Liu design
#'
#' Forward propagation over `(n, s)`.
#'
#' @param N Final sample size.
#' @param p True response rate.
#' @param table Output of [lee_liu_table()].
#' @param p0 Null response rate.
#' @param eta Final posterior probability threshold.
#' @param a,b Beta prior.
#' @return Named vector `c(P_succ, EN)`.
#' @export
#' @examples
#' tb <- lee_liu_table(20, 0.2, c(10, 15), 0.05, 0.95, 0.9)
#' lee_liu_oc(20, 0.2, tb, 0.2, 0.9)
lee_liu_oc <- function(N, p, table, p0, eta, a = 1, b = 1) {
  final_succ <- stats::pbeta(p0, a + 0:N, b + N - 0:N, lower.tail = FALSE) > eta
  cur <- 1
  P <- 0; EN <- 0
  for (n in seq_len(N)) {
    new <- numeric(n + 1)
    new[1:n] <- new[1:n] + cur * (1 - p)
    new[2:(n + 1)] <- new[2:(n + 1)] + cur * p
    cur <- new
    row <- which(table$n == n)
    if (length(row)) {
      sf <- table$s_fut[row]; se <- table$s_eff[row]
      for (s in 0:n) {
        if (cur[s + 1] == 0) next
        stop_f <- !is.na(sf) && s <= sf
        stop_e <- !is.na(se) && s >= se
        if (stop_e) {
          P <- P + cur[s + 1]; EN <- EN + n * cur[s + 1]; cur[s + 1] <- 0
        } else if (stop_f) {
          EN <- EN + n * cur[s + 1]; cur[s + 1] <- 0
        }
      }
    }
  }
  EN <- EN + N * sum(cur)
  P <- P + sum(cur[final_succ])
  c(P_succ = P, EN = EN)
}

#' Exact predictive probability, two arms
#'
#' Probability that `Pr(p1 > p0 | all data) > eta` after `m1`, `m0` further outcomes.
#' Posterior superiority is computed with the exact rank-scale integral of [prob_superior()].
#'
#' @param s1,n1,s0,n0 Responses and patients on treatment and control so far.
#' @param m1,m0 Further patients on each arm.
#' @param eta Threshold.
#' @param pri Beta priors `c(a1, b1, a0, b0)`.
#' @return A probability.
#' @export
#' @examples
#' pp_two_arm(8, 20, 5, 20, 10, 10, 0.9)
pp_two_arm <- function(s1, n1, s0, n0, m1, m0, eta, pri = c(1, 1, 1, 1)) {
  x1 <- 0:m1; x0 <- 0:m0
  w1 <- .dbetabinom(x1, m1, pri[1] + s1, pri[2] + n1 - s1)
  w0 <- .dbetabinom(x0, m0, pri[3] + s0, pri[4] + n0 - s0)
  Q <- .pr_superior_grid(pri[1] + s1 + x1, pri[2] + n1 - s1 + m1 - x1,
                         pri[3] + s0 + x0, pri[4] + n0 - s0 + m0 - x0)
  drop(w1 %*% (Q > eta) %*% w0)
}

#' Backward induction for a balanced two-arm binary trial
#'
#' Optimal stopping with per-pair sampling cost `c`. With `loss = "utility"` the terminal loss
#' is `-max(0, delta_hat)`; with `loss = "class"` it is `min(gamma (1 - q), (1 - gamma) q)` where
#' `q = Pr(p1 > p0 | data)`.
#'
#' @param T Maximum patients per arm.
#' @param c Cost per pair of patients.
#' @param loss `"utility"` or `"class"`.
#' @param gamma Classification-loss weight (also the declaration threshold on `q`).
#' @param pri Beta priors `c(a1, b1, a0, b0)`.
#' @return List with `stop` and `action`, each a list of length `T + 1` whose element `k + 1`
#'   is a `(k + 1) x (k + 1)` matrix indexed `[s1 + 1, s0 + 1]` (`action` 1 = declare treatment).
#' @export
#' @examples
#' bi <- backward_induction(10, 0.01); sum(bi$stop[[1]])
backward_induction <- function(T, c, loss = c("utility", "class"), gamma = 0.975,
                               pri = c(1, 1, 1, 1)) {
  loss <- match.arg(loss)
  a1 <- pri[1]; b1 <- pri[2]; a0 <- pri[3]; b0 <- pri[4]
  V <- vector("list", T + 1); stop <- V; act <- V
  for (k in T:0) {
    s <- 0:k
    A1 <- a1 + s; A0 <- a0 + s
    if (loss == "utility") {
      d <- outer(A1 / (a1 + b1 + k), A0 / (a0 + b0 + k), `-`)
      h <- -pmax(d, 0); a <- (d > 0) * 1L
    } else {
      q <- .pr_superior_grid(a1 + s, b1 + k - s, a0 + s, b0 + k - s)
      h <- pmin(gamma * (1 - q), (1 - gamma) * q); a <- (q > gamma) * 1L
    }
    if (k == T) {
      V[[k + 1]] <- h; stop[[k + 1]] <- matrix(TRUE, k + 1, k + 1); act[[k + 1]] <- a
      next
    }
    p1 <- matrix((a1 + s) / (a1 + b1 + k), k + 1, k + 1)
    p0 <- matrix((a0 + s) / (a0 + b0 + k), k + 1, k + 1, byrow = TRUE)
    Vn <- V[[k + 2]]
    i <- 1:(k + 1); j <- 1:(k + 1)
    cont <- c + ((1 - p1) * (1 - p0) * Vn[i, j] + p1 * (1 - p0) * Vn[i + 1, j] +
                 (1 - p1) * p0 * Vn[i, j + 1] + p1 * p0 * Vn[i + 1, j + 1])
    V[[k + 1]] <- pmin(h, cont); stop[[k + 1]] <- h <= cont; act[[k + 1]] <- a
  }
  list(stop = stop, action = act)
}

#' Exact operating characteristics of a backward-induction design
#'
#' Forward propagation of the count distribution.
#'
#' @param bi Output of [backward_induction()] (list with `stop`, `action`).
#' @param p1,p0 True response rates.
#' @return Named vector `c(P, EN)`: probability of declaring treatment, expected per-arm N.
#' @export
#' @examples
#' forward_oc(backward_induction(10, 0.01), 0.5, 0.3)
forward_oc <- function(bi, p1, p0) {
  stop <- bi$stop; act <- bi$action
  T <- length(stop) - 1
  dist <- matrix(1, 1, 1); P <- 0; EN <- 0
  for (k in 0:T) {
    st <- stop[[k + 1]]
    mass <- dist * st
    P <- P + sum(mass * act[[k + 1]]); EN <- EN + k * sum(mass)
    if (k == T) break
    dist <- dist * (!st)
    new <- matrix(0, k + 2, k + 2)
    i <- 1:(k + 1)
    new[i, i] <- new[i, i] + dist * (1 - p1) * (1 - p0)
    new[i + 1, i] <- new[i + 1, i] + dist * p1 * (1 - p0)
    new[i, i + 1] <- new[i, i + 1] + dist * (1 - p1) * p0
    new[i + 1, i + 1] <- new[i + 1, i + 1] + dist * p1 * p0
    dist <- new
  }
  c(P = P, EN = EN)
}

#' Exact operating characteristics of a fixed two-arm design
#'
#' Fixed design with `T` per arm declaring treatment if `Pr(p1 > p0 | data) > gamma`.
#'
#' @param T Per-arm sample size.
#' @param gamma Posterior probability threshold.
#' @param p1,p0 True response rates.
#' @param pri Beta priors `c(a1, b1, a0, b0)`.
#' @return Probability of declaring treatment.
#' @export
#' @examples
#' fixed_oc(20, 0.95, 0.5, 0.3)
fixed_oc <- function(T, gamma, p1, p0, pri = c(1, 1, 1, 1)) {
  s <- 0:T
  q <- .pr_superior_grid(pri[1] + s, pri[2] + T - s, pri[3] + s, pri[4] + T - s)
  sum(outer(stats::dbinom(s, T, p1), stats::dbinom(s, T, p0)) * (q > gamma))
}
