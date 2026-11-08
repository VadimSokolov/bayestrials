# Classical O'Brien-Fleming, 4 equally spaced looks, one-sided 0.025: gs_boundary(4)
.select_obf <- c(4.0487, 2.8628, 2.3375, 2.0243)
.select_fracs <- c(0.25, 0.5, 0.75, 1)

# SELECT-type reconstruction (Gaussian score process in information time).
# Information per event is 1/4, so I = d/4 and the score S has mean theta * I.

#' Posterior probability of benefit
#'
#' Pr(theta < 0 | data), theta = log HR, prior N(0, tau^2), score S, information I.
#'
#' @param S score statistic.
#' @param I information.
#' @param tau prior sd.
#' @return numeric.
#' @examples select_post_prob(-10, 238)
#' @export
select_post_prob <- function(S, I, tau = 0.25) {
  A <- I + 1 / tau^2
  pnorm(-S / sqrt(A))
}

#' Log Bayes factor (mixture e-value) for the score process
#'
#' @inheritParams select_post_prob
#' @return numeric.
#' @examples select_log_bf(-10, 238)
#' @export
select_log_bf <- function(S, I, tau = 0.25) {
  A <- I + 1 / tau^2
  -0.5 * log(1 + I * tau^2) + S^2 / (2 * A)
}

#' Monte Carlo of the calibrated SELECT rules
#'
#' Simulates Gaussian score paths at one-event steps of two, calibrates the posterior bar p* to type I error `alpha`, and compares O'Brien-Fleming (four looks), an e-value mixture and the Bayes rule. Paths are processed in batches to bound memory.
#'
#' @param d_tot total events.
#' @param HR true hazard ratio for the power calculation.
#' @param nsim number of simulated trials (per hypothesis).
#' @param tau prior sd for log HR.
#' @param alpha one-sided level.
#' @param seed optional seed.
#' @param batch trials per batch.
#' @return list with `naive_typeI`, `pstar`, `bayes_typeI`, `ev_typeI`, `gs_typeI`, `dg` (event grid) and `rows`, a data.frame with `rule`, `median_events`, `mean_events`, `power`.
#' @examples select_simulate(200, 0.8, nsim = 200, seed = 1)$pstar
#' @references Polson and Sokolov, case-study calibration for SELECT; Lincoff et al. (2023), NEJM 389, 2221-2232.
#' @export
select_simulate <- function(d_tot = 1270, HR = 0.80, nsim = 40000, tau = 0.25,
                            alpha = 0.025, seed = NULL, batch = 4000) {
  if (!is.null(seed)) set.seed(seed)
  B <- .select_obf
  fracs <- .select_fracs
  theta <- log(HR)
  dg <- seq(4, d_tot, by = 2)
  Ig <- dg / 4
  K <- length(dg)
  dI <- diff(c(0, Ig))
  kk <- vapply(fracs, function(f) which.min(abs(dg - f * d_tot)), 1L)
  sI <- sqrt(Ig)
  paths <- function(th, nb) {
    x <- matrix(rnorm(K * nb, sd = sqrt(dI)), K, nb)
    x <- apply(x, 2, cumsum)
    x + th * Ig
  }
  gs_hit <- function(S) {
    Z <- S[kk, , drop = FALSE] / sI[kk]
    colSums(Z <= -B) > 0
  }
  starts <- seq(1, nsim, by = batch)
  sizes <- pmin(batch, nsim - starts + 1)
  # pass 1: null
  naive <- ev0 <- gs0 <- logical(nsim)
  mx <- numeric(nsim)
  pos <- 0
  for (nb in sizes) {
    S <- paths(0, nb)
    pp <- select_post_prob(S, Ig, tau)
    idx <- pos + seq_len(nb)
    mx[idx] <- apply(pp, 2, max)
    ev0[idx] <- colSums(select_log_bf(S, Ig, tau) >= log(1 / alpha) & S < 0) > 0
    gs0[idx] <- gs_hit(S)
    pos <- pos + nb
  }
  kcut <- round(alpha * nsim)
  pstar <- sort(mx)[nsim - kcut + 1]
  # pass 2: alternative
  first <- list(gs = rep(NA_real_, nsim), ev = rep(NA_real_, nsim), bayes = rep(NA_real_, nsim))
  pos <- 0
  firstidx <- function(m) {
    hit <- colSums(m) > 0
    ifelse(hit, dg[apply(m, 2, which.max)], NA_real_)
  }
  for (nb in sizes) {
    S <- paths(theta, nb)
    idx <- pos + seq_len(nb)
    Z <- S / sI
    g <- matrix(FALSE, K, nb)
    for (j in seq_along(kk)) g[kk[j], ] <- Z[kk[j], ] <= -B[j]
    first$gs[idx] <- firstidx(g)
    first$ev[idx] <- firstidx(select_log_bf(S, Ig, tau) >= log(1 / alpha) & S < 0)
    first$bayes[idx] <- firstidx(select_post_prob(S, Ig, tau) >= pstar)
    pos <- pos + nb
  }
  nm <- c("O'Brien-Fleming, four scheduled looks",
          "E-value mixture, continuous, E >= 40",
          sprintf("Bayes posterior, continuous, p* = %.3f", pstar))
  rows <- data.frame(
    rule = nm,
    median_events = vapply(first, median, 0, na.rm = TRUE),
    mean_events = vapply(first, mean, 0, na.rm = TRUE),
    power = vapply(first, function(x) mean(!is.na(x)), 0),
    row.names = NULL)
  # pass-1 naive error at the nominal 0.975 posterior bar
  naive_typeI <- mean(mx >= 0.975)
  list(naive_typeI = naive_typeI, pstar = pstar, bayes_typeI = mean(mx >= pstar),
       ev_typeI = mean(ev0), gs_typeI = mean(gs0), dg = dg, rows = rows)
}

#' Calendar month at which a fraction of events has accrued
#'
#' Linear enrolment over `T_en` months, events proportional to person time.
#'
#' @param d events at the decision.
#' @param d_tot total events.
#' @param T_en enrolment duration in months.
#' @param T_end trial duration in months.
#' @return month (numeric, vectorized over `d`).
#' @examples select_month_of_events(464, 1270)
#' @export
select_month_of_events <- function(d, d_tot, T_en = 36, T_end = 60) {
  pt <- function(t) ifelse(t <= T_en, 0.5 * t^2 / T_en, 0.5 * T_en + (t - T_en))
  one <- function(d) uniroot(function(t) pt(t) / pt(T_end) - d / d_tot,
                             c(1e-9, T_end), tol = 1e-10)$root
  vapply(d, one, 0)
}

#' Look table at the observed effect
#'
#' @param d_tot total events.
#' @param HR hazard ratio assumed observed at each look.
#' @param tau prior sd.
#' @return data.frame with `events`, `z`, `obf_bound`, `log_bf`, `post_prob`.
#' @examples select_look_table(952, 0.8)
#' @export
select_look_table <- function(d_tot, HR, tau = 0.25) {
  f <- .select_fracs
  B <- .select_obf
  theta <- log(HR)
  d <- f * d_tot
  I <- d / 4
  S <- theta * I
  data.frame(events = as.integer(round(d)), z = abs(S) / sqrt(I), obf_bound = B,
             log_bf = select_log_bf(S, I, tau), post_prob = select_post_prob(S, I, tau))
}
