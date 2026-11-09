#' Example dose-toxicity scenarios
#'
#' @return A named list of true DLT probability vectors (A, B, C) for six doses.
#' @examples dose_scenarios()$A
#' @export
dose_scenarios <- function() {
  list(A = c(0.05, 0.10, 0.25, 0.40, 0.55, 0.70),
       B = c(0.02, 0.05, 0.08, 0.12, 0.25, 0.45),
       C = c(0.10, 0.25, 0.45, 0.60, 0.70, 0.80))
}

#' True maximum tolerated dose
#'
#' @param p True DLT probabilities by dose.
#' @param target Target toxicity rate.
#' @return Highest dose (1-based) with rate at most the target, 0 if none.
#' @examples dose_true_mtd(c(0.05, 0.1, 0.25, 0.4))
#' @export
dose_true_mtd <- function(p, target = 0.25) {
  ok <- which(p <= target + 1e-12)
  if (length(ok)) max(ok) else 0L
}

#' One simulated 3+3 trial
#'
#' The rule implemented, step by step. Treat a cohort of 3 at dose 1. At the current dose: with 0 DLTs
#' in 3, escalate; with 1 DLT in 3, treat 3 more at the same dose, escalate if the total is 1 DLT in 6
#' and otherwise stop; with 2 or 3 DLTs in the first 3, stop immediately (no expansion cohort). Escalation
#' moves up exactly one dose, and a dose is never revisited: there is no de-escalation to a lower dose
#' for further testing, no expansion of the dose below the stopping dose, and no confirmation of the
#' selected dose. When the trial stops for excess toxicity at dose `d` the recommended dose is `d - 1`
#' (0, "none", if `d = 1`), which may have been given to only 3 patients. If the top dose passes (0 DLTs
#' in 3, or 1 DLT in 6 after expansion) the top dose is selected, with no expansion cohort after 0 DLTs
#' in 3. There is no sample size cap, so the trial length is random. This is the traditional "3+3" of the usual
#' textbook form; Storer's (1989) designs, which expand and de-escalate after initial toxicity, are
#' different rules.
#'
#' @param p True DLT probabilities by dose.
#' @return List with `selected` (recommended dose, 1-based, 0 if none), `n` and `y` (patients and DLTs per dose).
#' @references Storer, B. E. (1989). Design and analysis of phase I clinical trials. Biometrics 45, 925-937.
#' @examples set.seed(1); dose_three_plus_three(c(0.05, 0.1, 0.25, 0.4))$selected
#' @export
dose_three_plus_three <- function(p) {
  k <- length(p)
  n <- integer(k); y <- integer(k)
  d <- 1L
  repeat {
    a <- stats::rbinom(1, 3, p[d]); n[d] <- n[d] + 3L; y[d] <- y[d] + a
    if (a == 0) {
      if (d == k) return(list(selected = k, n = n, y = y))
      d <- d + 1L; next
    }
    if (a == 1) {
      b <- stats::rbinom(1, 3, p[d]); n[d] <- n[d] + 3L; y[d] <- y[d] + b
      if (a + b == 1) {
        if (d == k) return(list(selected = k, n = n, y = y))
        d <- d + 1L; next
      }
    }
    return(list(selected = d - 1L, n = n, y = y))
  }
}

#' One simulated one-parameter power-model CRM trial
#'
#' Dose-toxicity model `p_d(a) = skeleton_d^exp(a)` with a N(0, 1.34) prior on `a`. The posterior is
#' computed on a grid of 401 equally spaced values of `a` on [-4, 4] (the prior and posterior are
#' truncated to this range and renormalised on it). The trial starts at dose 1 and enrols `N %/% cohort`
#' cohorts of size `cohort` (patients beyond `cohort * (N %/% cohort)` are not enrolled). After each
#' cohort the posterior is updated and: (i) the trial stops with no dose selected (`selected = 0`) if
#' the posterior probability that the lowest dose's toxicity exceeds `target` is more than 0.90 (this
#' safety cutoff is checked after every cohort, including the first, with no minimum number of
#' patients); otherwise (ii) the recommended dose is the one whose posterior mean toxicity is closest to
#' `target` (either side of the target, so the recommended dose can exceed the target) and (iii) the next
#' cohort is treated at that dose, except that escalation never skips: the next dose is at most one above
#' the current one, while de-escalation can jump several levels. The selected dose is the recommendation
#' computed after the last cohort. It is not restricted by the no-skipping rule and may be a dose that was
#' never given. There is no early stopping for sample size or for a precisely estimated MTD.
#'
#' @param p True DLT probabilities by dose.
#' @param skeleton Prior guesses of DLT probabilities, one per dose (same length as `p`).
#' @param N Maximum sample size.
#' @param cohort Cohort size.
#' @param target Target toxicity rate.
#' @return List with `selected`, `n`, `y` as in [dose_three_plus_three()].
#' @references O'Quigley, J., Pepe, M. and Fisher, L. (1990). Continual reassessment method. Biometrics 46, 33-48.
#' @examples set.seed(1); dose_crm(c(0.05, 0.1, 0.25, 0.4, 0.55, 0.7), N = 9)$selected
#' @export
dose_crm <- function(p, skeleton = c(0.05, 0.12, 0.25, 0.40, 0.55, 0.70), N = 30, cohort = 3,
                     target = 0.25) {
  k <- length(p)
  if (length(skeleton) != k) stop("`skeleton` must have the same length as `p`")
  grid <- seq(-4, 4, length.out = 401)
  prior <- stats::dnorm(grid, 0, sqrt(1.34)); prior <- prior / sum(prior)
  pw <- outer(skeleton, exp(grid), `^`)
  n <- integer(k); y <- integer(k)
  d <- 1L; new <- 1L
  for (i in seq_len(N %/% cohort)) {
    n[d] <- n[d] + as.integer(cohort); y[d] <- y[d] + stats::rbinom(1, cohort, p[d])
    ll <- colSums(y * log(pw) + (n - y) * log1p(-pw))
    post <- prior * exp(ll - max(ll)); post <- post / sum(post)
    if (sum(post * (pw[1, ] > target)) > 0.90) return(list(selected = 0L, n = n, y = y))
    est <- as.vector(pw %*% post)
    new <- which.min(abs(est - target))
    d <- min(new, d + 1L)
  }
  list(selected = new, n = n, y = y)
}

#' BOIN escalation and de-escalation boundaries
#'
#' @param target Target toxicity rate.
#' @param f1,f2 Multipliers defining the sub-therapeutic (f1 * target) and overly toxic (f2 * target) rates.
#' @return Named numeric vector with `escalate` and `deescalate` cutoffs for the observed DLT rate.
#' @references Liu, S. and Yuan, Y. (2015). Bayesian optimal interval designs for phase I clinical trials. Journal of the Royal Statistical Society C 64, 507-523.
#' @examples dose_boin_bounds(0.25)
#' @export
dose_boin_bounds <- function(target = 0.25, f1 = 0.6, f2 = 1.4) {
  p1 <- f1 * target; p2 <- f2 * target
  le <- log((1 - p1) / (1 - target)) / log(target * (1 - p1) / (p1 * (1 - target)))
  ld <- log((1 - target) / (1 - p2)) / log(p2 * (1 - target) / (target * (1 - p2)))
  c(escalate = le, deescalate = ld)
}

#' Pool adjacent violators (weighted isotonic regression, increasing)
#'
#' @param r Raw rates.
#' @param w Weights.
#' @return Isotonic fit, same length as `r`.
#' @examples dose_pava(c(0.1, 0.3, 0.2, 0.5), c(1, 1, 2, 1))
#' @export
dose_pava <- function(r, w) {
  val <- numeric(0); wt <- numeric(0); cnt <- integer(0); m <- 0L
  for (i in seq_along(r)) {
    m <- m + 1L; val[m] <- r[i]; wt[m] <- w[i]; cnt[m] <- 1L
    while (m > 1L && val[m - 1L] > val[m]) {
      tw <- wt[m - 1L] + wt[m]
      val[m - 1L] <- (val[m - 1L] * wt[m - 1L] + val[m] * wt[m]) / tw
      wt[m - 1L] <- tw; cnt[m - 1L] <- cnt[m - 1L] + cnt[m]
      m <- m - 1L
    }
  }
  rep(val[seq_len(m)], cnt[seq_len(m)])
}

# Final BOIN dose selection, following BOIN::select.mtd (Liu and Yuan 2015, as implemented by the
# authors' package): doses are eliminated from the first dose i with n_i >= 3 and
# P(p_i > target | data) > cutoff (Beta(1 + y, 1 + n - y) posterior) upwards, recomputed from the
# final data; among the tried, non-eliminated doses the isotonic regression (PAVA) of
# (y + 0.05) / (n + 0.1) is fitted with weights equal to the inverse posterior variance
# 1 / [(y + 0.05)(n - y + 0.05) / ((n + 0.1)^2 (n + 1.1))]; the dose with fitted rate closest to the
# target is selected. Ties are broken as in the reference, which adds i * 1e-10 to the i-th fitted rate:
# the higher dose among tied doses below the target, the lower dose among tied doses above it.
# Returns 0 if the first dose is eliminated or no non-eliminated dose has data.
dose_boin_select <- function(n, y, target = 0.25, cutoff = 0.95) {
  k <- length(n)
  elim <- rep(FALSE, k)
  for (i in seq_len(k)) {
    if (n[i] >= 3 && stats::pbeta(target, y[i] + 1, n[i] - y[i] + 1, lower.tail = FALSE) > cutoff) {
      elim[i:k] <- TRUE; break
    }
  }
  if (elim[1] || sum(n[!elim]) == 0) return(0L)
  adm <- which(n != 0 & !elim)
  ya <- y[adm]; na <- n[adm]
  phat <- (ya + 0.05) / (na + 0.1)
  v <- (ya + 0.05) * (na - ya + 0.05) / ((na + 0.1)^2 * (na + 1.1))
  iso <- dose_pava(phat, 1 / v) + seq_along(adm) * 1e-10
  as.integer(adm[order(abs(iso - target))[1]])
}

#' One simulated BOIN trial
#'
#' Dosing follows the BOIN interval rule: with the observed DLT rate `y/n` at the current dose, escalate
#' when it is at most the escalation boundary, de-escalate when it is at least the de-escalation boundary
#' (see [dose_boin_bounds()]) and otherwise stay; escalation never goes to or beyond an eliminated dose.
#' A dose and all higher doses are eliminated once `n >= 3` and `P(p > target | data) > 0.95` under a
#' Beta(1, 1) prior (elimination is permanent during the trial); if dose 1 is eliminated the trial stops
#' with `selected = 0`. There is no early stop for large `n` at a dose and no safety cap other than
#' elimination. Patients beyond `cohort * (N %/% cohort)` are not enrolled.
#'
#' Final selection reproduces the reference `select.mtd` of the BOIN package: the elimination rule is
#' re-evaluated on the final data, and isotonic regression with inverse posterior-variance weights
#' (not patient counts) is applied to `(y + 0.05) / (n + 0.1)` over the tried, non-eliminated doses. The
#' dose with fitted rate closest to the target is selected. Ties are broken as in the reference: the
#' highest dose among tied fitted rates below the target and the lowest dose among tied fitted rates
#' above it. (Earlier versions of this function used patient-count weights and always picked the lowest
#' tied dose.)
#'
#' @param p True DLT probabilities by dose.
#' @param N Maximum sample size.
#' @param cohort Cohort size.
#' @param target Target toxicity rate.
#' @return List with `selected`, `n`, `y` as in [dose_three_plus_three()].
#' @references Liu, S. and Yuan, Y. (2015). Bayesian optimal interval designs for phase I clinical trials. Journal of the Royal Statistical Society C 64, 507-523.
#' @examples set.seed(1); dose_boin(c(0.05, 0.1, 0.25, 0.4, 0.55, 0.7), N = 9)$selected
#' @export
dose_boin <- function(p, N = 30, cohort = 3, target = 0.25) {
  k <- length(p)
  b <- dose_boin_bounds(target); le <- b[["escalate"]]; ld <- b[["deescalate"]]
  n <- integer(k); y <- integer(k)
  d <- 1L; top <- k + 1L  # doses with index >= top are eliminated
  for (i in seq_len(N %/% cohort)) {
    n[d] <- n[d] + as.integer(cohort); y[d] <- y[d] + stats::rbinom(1, cohort, p[d])
    for (j in seq_len(k)) {
      if (n[j] >= 3 && stats::pbeta(target, 1 + y[j], 1 + n[j] - y[j], lower.tail = FALSE) > 0.95) {
        top <- min(top, j); break
      }
    }
    if (top == 1L) return(list(selected = 0L, n = n, y = y))
    r <- y[d] / n[d]
    if (r <= le) d <- min(d + 1L, top - 1L)
    else if (r >= ld) d <- max(d - 1L, 1L)
  }
  list(selected = dose_boin_select(n, y, target), n = n, y = y)
}

#' Operating characteristics of a dose-escalation design by simulation
#'
#' @param design Function of the true DLT probability vector returning a list with `selected`, `n`, `y`
#'   (e.g. [dose_three_plus_three], [dose_crm], [dose_boin]).
#' @param p True DLT probabilities.
#' @param reps Number of simulated trials.
#' @param seed Optional seed.
#' @param target Target toxicity rate defining the true MTD.
#' @param toxic Doses with true rate above this count as overly toxic.
#' @return List: `p_correct`, `p_over`, `p_under`, mean `n`, `n_toxic`, `dlts`, and selection probabilities `sel` (element 1 is "none", then doses 1..K).
#' @examples dose_simulate(dose_three_plus_three, c(0.05, 0.1, 0.25, 0.4), reps = 50, seed = 1)$p_correct
#' @export
dose_simulate <- function(design, p, reps, seed = NULL, target = 0.25, toxic = 0.35) {
  if (!is.null(seed)) set.seed(seed)
  m <- dose_true_mtd(p, target)
  sel <- numeric(length(p) + 1)
  nt <- dl <- ntox <- numeric(reps)
  for (i in seq_len(reps)) {
    r <- design(p)
    sel[r$selected + 1] <- sel[r$selected + 1] + 1
    nt[i] <- sum(r$n); dl[i] <- sum(r$y); ntox[i] <- sum(r$n[p > toxic])
  }
  sel <- sel / reps
  list(p_correct = sel[m + 1], p_over = sum(sel[-seq_len(m + 1)]),
       p_under = sum(sel[seq_len(m)]), n = mean(nt), n_toxic = mean(ntox),
       dlts = mean(dl), sel = sel)
}
