#' Simulate one confounded trial
#'
#' Hidden health `h ~ N(0,1)` raises recovery. With `assign = "doctor"` the drug goes preferentially to healthier patients; with a coin it ignores `h`.
#'
#' @param n Number of patients.
#' @param effect True additive effect of the drug on recovery probability.
#' @param assign `"doctor"` or `"coin"`.
#' @param seed Optional seed.
#' @return Named vector: `naive` (naive difference in recovery) and `health_gap` (mean health, drug minus control).
#' @examples
#' sim_confounded_trial(200, seed = 1)
#' @export
sim_confounded_trial <- function(n, effect = 0.10, assign = c("doctor", "coin"), seed = NULL) {
  assign <- match.arg(assign)
  if (!is.null(seed)) set.seed(seed)
  h <- stats::rnorm(n)
  p_drug <- if (assign == "doctor") stats::plogis(1.5 * h) else rep(0.5, n)
  d <- stats::runif(n) < p_drug
  p0 <- stats::plogis(-0.2 + 0.9 * h)
  p1 <- pmin(pmax(p0 + effect, 0), 1)
  y <- stats::runif(n) < ifelse(d, p1, p0)
  c(naive = mean(y[d]) - mean(y[!d]), health_gap = mean(h[d]) - mean(h[!d]))
}

#' Doctor-chosen trial analysed naively and stratified on health bands
#'
#' @param n Number of patients.
#' @param effect True additive effect of the drug.
#' @param bands Number of quantile bands of the recorded health score.
#' @param seed Optional seed.
#' @return Named vector: `adjusted` (band-stratified, size-weighted) and `naive`.
#' @examples
#' sim_adjusted_vs_naive(500, seed = 1)
#' @export
sim_adjusted_vs_naive <- function(n, effect = 0.10, bands = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  h <- stats::rnorm(n)
  d <- stats::runif(n) < stats::plogis(1.5 * h)
  p0 <- stats::plogis(-0.2 + 0.9 * h)
  y <- stats::runif(n) < ifelse(d, pmin(pmax(p0 + effect, 0), 1), p0)
  q <- stats::quantile(h, seq(0, 1, length.out = bands + 1), names = FALSE)
  g <- pmin(pmax(findInterval(h, q, left.open = FALSE) , 1), bands)
  est <- numeric(0); w <- numeric(0)
  for (k in seq_len(bands)) {
    m <- g == k
    if (any(d & m) && any(!d & m)) {
      est <- c(est, mean(y[d & m]) - mean(y[!d & m]))
      w <- c(w, sum(m))
    }
  }
  c(adjusted = stats::weighted.mean(est, w), naive = mean(y[d]) - mean(y[!d]))
}

#' Simple (coin-flip) assignment
#'
#' @param n Number of patients.
#' @param seed Optional seed.
#' @return Integer vector of 0/1 arms.
#' @examples
#' assign_simple(10, seed = 1)
#' @export
assign_simple <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  sample.int(2, n, replace = TRUE) - 1L
}

#' Permuted-block assignment
#'
#' @param n Number of patients.
#' @param sizes Even block sizes, one drawn at random per block.
#' @param seed Optional seed.
#' @return Integer vector of 0/1 arms.
#' @examples
#' assign_block(12, seed = 1)
#' @export
assign_block <- function(n, sizes = 4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- integer(0)
  while (length(out) < n) {
    b <- sizes[sample.int(length(sizes), 1)]
    blk <- rep(c(0L, 1L), b %/% 2)
    out <- c(out, blk[sample.int(length(blk))])
  }
  out[seq_len(n)]
}

#' Stratified permuted-block assignment
#'
#' Permuted blocks run separately inside each stratum.
#'
#' @param strata Vector of stratum labels.
#' @param sizes Even block sizes.
#' @param seed Optional seed.
#' @return Integer vector of 0/1 arms.
#' @examples
#' assign_stratified(rep(1:2, 6), seed = 1)
#' @export
assign_stratified <- function(strata, sizes = 4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  out <- integer(length(strata))
  for (s in sort(unique(strata))) {
    idx <- which(strata == s)
    out[idx] <- assign_block(length(idx), sizes)
  }
  out
}

#' Pocock-Simon minimization
#'
#' Marginal minimization over the columns of `X`: assign to the arm with the smaller count of same-level patients summed over factors, with probability `p_follow`.
#'
#' @param X Matrix of integer factor levels (rows are patients in order of arrival).
#' @param p_follow Probability of following the minimization rule.
#' @param seed Optional seed.
#' @return Integer vector of 0/1 arms.
#' @references Pocock SJ, Simon R (1975) Biometrics 31:103-115.
#' @examples
#' assign_minimization(cbind(rep(0:1, 10), rep(0:1, each = 10)), seed = 1)
#' @export
assign_minimization <- function(X, p_follow = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X); m <- ncol(X)
  arm <- integer(n)
  counts <- vector("list", 2)
  counts[[1]] <- new.env(); counts[[2]] <- new.env()
  key <- function(j, v) paste(j, v, sep = "_")
  get_ct <- function(k, kk) { v <- counts[[k]][[kk]]; if (is.null(v)) 0 else v }
  for (i in seq_len(n)) {
    ks <- vapply(seq_len(m), function(j) key(j, X[i, j]), character(1))
    sc <- c(sum(vapply(ks, function(kk) get_ct(1, kk), numeric(1))),
            sum(vapply(ks, function(kk) get_ct(2, kk), numeric(1))))
    if (sc[1] == sc[2]) {
      a <- sample.int(2, 1) - 1L
    } else {
      best <- which.min(sc) - 1L
      a <- if (stats::runif(1) < p_follow) best else 1L - best
    }
    arm[i] <- a
    for (kk in ks) counts[[a + 1]][[kk]] <- get_ct(a + 1, kk) + 1
  }
  arm
}

#' Covariate imbalance between arms
#'
#' @param arm 0/1 arm vector.
#' @param x 0/1 covariate.
#' @return Absolute difference in proportion with x = 1 between arms.
#' @examples
#' covariate_imbalance(c(0, 1, 0, 1), c(1, 1, 0, 0))
#' @export
covariate_imbalance <- function(arm, x) abs(mean(x[arm == 1]) - mean(x[arm == 0]))

#' One-sided permutation p-value given an assignment
#'
#' Permutes `arm` (optionally within restricted blocks) and compares the difference in means.
#'
#' @param y Outcomes.
#' @param arm 0/1 arm vector.
#' @param nsim Number of permutations.
#' @param restricted_blocks Optional list of index vectors; permute within each.
#' @param seed Optional seed.
#' @return Named vector with `obs` and `p`.
#' @examples
#' perm_test_arm(c(1, 2, 3, 4, 5, 6), c(0, 0, 0, 1, 1, 1), nsim = 200, seed = 1)
#' @export
perm_test_arm <- function(y, arm, nsim = 5000, restricted_blocks = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  obs <- mean(y[arm == 1]) - mean(y[arm == 0])
  cnt <- 0
  for (i in seq_len(nsim)) {
    if (is.null(restricted_blocks)) {
      a <- arm[sample.int(length(arm))]
    } else {
      a <- arm
      for (blk in restricted_blocks) a[blk] <- arm[blk][sample.int(length(blk))]
    }
    cnt <- cnt + ((mean(y[a == 1]) - mean(y[a == 0])) >= obs)
  }
  c(obs = obs, p = (cnt + 1) / (nsim + 1))
}

#' Arnold stochasticity statistic
#'
#' @param points Event times on a circle of circumference `L`.
#' @param L Circumference.
#' @return Named vector: `S` (T times sum of squared normalised gaps) and `expected` (2T/(T+1)).
#' @examples
#' arnold_stochasticity(c(1, 3, 6), 10)
#' @export
arnold_stochasticity <- function(points, L) {
  pts <- sort(points %% L)
  T <- length(pts)
  d <- diff(c(pts, pts[1] + L))
  c(S = T * sum((d / L)^2), expected = 2 * T / (T + 1))
}

#' Compression ratio of a bit sequence
#'
#' Ratio of zlib-compressed to raw length of the sequence written as a string of digits. Uses `memCompress(type = "gzip")`, which is zlib at the default level (the Python reference uses level 9), so values can differ slightly.
#'
#' @param bits Vector of 0/1.
#' @return Numeric compression ratio.
#' @examples
#' compress_ratio(rep(0:1, 50))
#' @export
compress_ratio <- function(bits) {
  s <- charToRaw(paste(bits, collapse = ""))
  length(memCompress(s, type = "gzip")) / length(s)
}

#' Bits from a multiplicative congruential generator
#'
#' Iterates x_{k+1} = a x_k mod m from x_1 = `seed` and returns the indicator
#' x_k > m/2 rounded down. The default is the 16-bit generator 75 x mod 65537
#' (the ZX81 generator), a short formula whose bits look random to a compressor.
#'
#' @param n Number of bits.
#' @param seed Starting state x_1.
#' @param a Multiplier.
#' @param m Modulus.
#' @return An integer vector of 0/1 of length `n`.
#' @examples rand_bits_lcg(10)
#' @export
rand_bits_lcg <- function(n, seed = 1, a = 75, m = 65537) {
  x <- numeric(n); x[1] <- seed
  for (i in seq_len(n)[-1]) x[i] <- (a * x[i - 1]) %% m
  as.integer(x > floor(m / 2))
}
