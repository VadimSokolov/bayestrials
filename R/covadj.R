#' Variance and sample-size gain from ANCOVA adjustment
#'
#' Adjusting for a baseline covariate that explains a fraction `r2` of the outcome variance multiplies the variance of the estimated treatment effect by about `1 - r2`, and the sample size needed for the same power by the same factor.
#'
#' @param r2 Squared correlation (R-squared) between the baseline covariate(s) and the outcome within arms, in `[0, 1)`.
#' @param n_unadjusted Optional sample size per arm of the unadjusted analysis.
#' @return A list with `var_factor`, `n_ratio` and, if `n_unadjusted` is given, `n_adjusted` (rounded up).
#' @examples
#' covadj_ancova_gain(0.30, n_unadjusted = 133)
#' @export
covadj_ancova_gain <- function(r2, n_unadjusted = NULL) {
  stopifnot(is.numeric(r2), all(r2 >= 0), all(r2 < 1))
  out <- list(var_factor = 1 - r2, n_ratio = 1 - r2)
  if (!is.null(n_unadjusted)) out$n_adjusted <- ceiling((1 - r2) * n_unadjusted)
  out
}

#' Variance of the treatment effect under three analyses of a baseline-and-follow-up outcome
#'
#' With baseline-to-follow-up correlation `rho` and equal variances, the variance of the estimated treatment effect relative to the final-value analysis is 1 for final values alone, `2 (1 - rho)` for the change score and `1 - rho^2` for ANCOVA.
#'
#' @param rho Correlation between baseline and follow-up measurement.
#' @return Data frame with columns `rho`, `final`, `change`, `ancova`.
#' @examples
#' covadj_variance_factors(c(0.3, 0.5, 0.7))
#' @export
covadj_variance_factors <- function(rho) {
  stopifnot(is.numeric(rho), all(abs(rho) < 1))
  data.frame(rho = rho, final = 1, change = 2 * (1 - rho), ancova = 1 - rho^2)
}

#' True marginal log odds ratio implied by a logistic model with one normal covariate
#'
#' For `logit P(Y = 1) = a + b x + g trt` with `x ~ N(0, 1)`, the marginal (population-averaged) log odds ratio is smaller in absolute value than the conditional `g` when `b != 0`: the odds ratio is non-collapsible.
#'
#' @param a Intercept.
#' @param b Coefficient of the covariate.
#' @param g Conditional log odds ratio of treatment.
#' @return Named vector with the marginal risks `p1`, `p0`, the risk difference `rd` and the marginal log odds ratio `log_or`.
#' @examples
#' covadj_true_marginal(-0.5, 1.5, 0.7)
#' @export
covadj_true_marginal <- function(a, b, g) {
  m <- function(eta0) stats::integrate(function(x) stats::plogis(eta0 + b * x) * stats::dnorm(x), -Inf, Inf)$value
  p0 <- m(a); p1 <- m(a + g)
  c(p1 = p1, p0 = p0, rd = p1 - p0, log_or = stats::qlogis(p1) - stats::qlogis(p0))
}

#' Simulate a randomized trial with a prognostic baseline covariate and a binary outcome
#'
#' Treatment is assigned by complete randomization with equal arms, so it is independent of the covariate by design.
#'
#' @param n Total sample size (rounded down to an even number).
#' @param a,b,g Intercept, covariate coefficient and conditional log odds ratio of `logit P(Y = 1) = a + b x + g trt`.
#' @param seed Optional seed.
#' @return Data frame with `trt` (0/1), `x` (standard normal) and `y` (0/1).
#' @examples
#' head(covadj_sim_trial(100, seed = 1))
#' @export
covadj_sim_trial <- function(n, a = -0.5, b = 1.5, g = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- 2 * (n %/% 2)
  trt <- sample(rep(0:1, each = n / 2))
  x <- stats::rnorm(n)
  y <- stats::rbinom(n, 1, stats::plogis(a + b * x + g * trt))
  data.frame(trt = trt, x = x, y = y)
}

#' Standardized (g-computation) marginal effects from a fitted logistic model
#'
#' Predicts every patient's risk under treatment and under control from the fitted model, averages over the trial sample, and reports the marginal risk difference, relative risk and odds ratio. The standard error uses the influence function of the plug-in estimator, which accounts for the sampling variation of the covariates as well as of the coefficients (model-based coefficient covariance).
#'
#' @param fit A `glm` fit with `family = binomial()` (logit link) containing the treatment indicator.
#' @param data The data frame the model was fitted to.
#' @param trt Name of the 0/1 treatment column.
#' @return Data frame with one row per measure (`rd`, `rr`, `or`) and columns `estimate`, `se` (for `rr` and `or` on the log scale), `lower`, `upper` (95 percent interval on the natural scale). Attribute `risks` holds the marginal risks under each arm.
#' @examples
#' d <- covadj_sim_trial(400, seed = 2)
#' covadj_standardize(stats::glm(y ~ trt + x, binomial, d), d, "trt")
#' @export
covadj_standardize <- function(fit, data, trt) {
  stopifnot(inherits(fit, "glm"), fit$family$family == "binomial", fit$family$link == "logit")
  n <- nrow(data)
  d1 <- d0 <- data; d1[[trt]] <- 1; d0[[trt]] <- 0
  X1 <- stats::model.matrix(stats::formula(fit), d1)
  X0 <- stats::model.matrix(stats::formula(fit), d0)
  X <- stats::model.matrix(fit)
  b <- stats::coef(fit)
  p1 <- stats::plogis(drop(X1 %*% b)); p0 <- stats::plogis(drop(X0 %*% b))
  m1 <- mean(p1); m0 <- mean(p0)
  beta_if <- n * (X * (fit$y - stats::fitted(fit))) %*% stats::vcov(fit)
  if1 <- (p1 - m1) + drop(beta_if %*% colMeans(p1 * (1 - p1) * X1))
  if0 <- (p0 - m0) + drop(beta_if %*% colMeans(p0 * (1 - p0) * X0))
  se_of <- function(g1, g0) sqrt(mean((g1 * if1 + g0 * if0)^2) / n)
  est <- c(rd = m1 - m0, rr = log(m1 / m0), or = stats::qlogis(m1) - stats::qlogis(m0))
  se <- c(rd = se_of(1, -1), rr = se_of(1 / m1, -1 / m0),
          or = se_of(1 / (m1 * (1 - m1)), -1 / (m0 * (1 - m0))))
  z <- stats::qnorm(0.975)
  lo <- est - z * se; hi <- est + z * se
  tr <- c(rd = identity, rr = exp, or = exp)
  out <- data.frame(measure = names(est),
                    estimate = mapply(function(f, v) f(v), tr, est),
                    se = se,
                    lower = mapply(function(f, v) f(v), tr, lo),
                    upper = mapply(function(f, v) f(v), tr, hi),
                    row.names = NULL)
  attr(out, "risks") <- c(treated = m1, control = m0)
  out
}

#' Repeated simulated trials: unadjusted, covariate-adjusted conditional and standardized marginal log odds ratios
#'
#' For each simulated trial fits the unadjusted logistic model, the adjusted logistic model (conditional log odds ratio) and standardizes the adjusted fit (marginal log odds ratio).
#'
#' @param nsim Number of trials.
#' @param n Total sample size per trial.
#' @param a,b,g Parameters of `covadj_sim_trial()`.
#' @param seed Seed.
#' @return Data frame with, for each of `unadj`, `adj`, `std`, the estimate (`_est`) and standard error (`_se`) of the log odds ratio.
#' @examples
#' head(covadj_sim_compare(20, 200, seed = 1))
#' @export
covadj_sim_compare <- function(nsim, n, a = -0.5, b = 1.5, g = 0.7, seed = 1) {
  set.seed(seed)
  res <- t(vapply(seq_len(nsim), function(i) {
    d <- covadj_sim_trial(n, a, b, g)
    f0 <- suppressWarnings(stats::glm(y ~ trt, stats::binomial(), d))
    f1 <- suppressWarnings(stats::glm(y ~ trt + x, stats::binomial(), d))
    s <- covadj_standardize(f1, d, "trt")
    s <- s[s$measure == "or", ]
    c(stats::coef(f0)[["trt"]], sqrt(stats::vcov(f0)["trt", "trt"]),
      stats::coef(f1)[["trt"]], sqrt(stats::vcov(f1)["trt", "trt"]),
      log(s$estimate), s$se)
  }, numeric(6)))
  colnames(res) <- c("unadj_est", "unadj_se", "adj_est", "adj_se", "std_est", "std_se")
  as.data.frame(res)
}

#' Bayesian g-computation for a binary outcome: posterior of the marginal effect
#'
#' Fits a logistic model with independent normal priors on the coefficients (centered at zero, standard deviation `prior_sd`; a wide `N(0, 10^2)` prior on the intercept), approximates the posterior by a normal distribution at the mode (Laplace), and for each posterior draw averages the predicted risks under treatment and under control over the trial sample. The draws of the averages are draws from the posterior of the marginal risk difference and log odds ratio.
#'
#' @param formula Model formula with the treatment indicator and any baseline covariates.
#' @param data Data frame.
#' @param trt Name of the 0/1 treatment column.
#' @param prior_sd Prior standard deviation of the non-intercept coefficients.
#' @param ndraw Number of posterior draws.
#' @param seed Optional seed.
#' @return Data frame of draws with columns `m1`, `m0`, `rd` and `log_or`.
#' @examples
#' d <- covadj_sim_trial(200, seed = 3)
#' dr <- covadj_bayes_gcomp(y ~ trt + x, d, "trt", seed = 1)
#' c(mean = mean(dr$rd), sd = stats::sd(dr$rd))
#' @export
covadj_bayes_gcomp <- function(formula, data, trt, prior_sd = 2.5, ndraw = 4000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- stats::model.matrix(formula, data)
  y <- stats::model.response(stats::model.frame(formula, data))
  psd <- c(10, rep(prior_sd, ncol(X) - 1))
  nlp <- function(b) {
    eta <- drop(X %*% b)
    -sum(y * eta - log1p(exp(eta))) + sum(b^2 / (2 * psd^2))
  }
  gr <- function(b) {
    drop(crossprod(X, stats::plogis(drop(X %*% b)) - y)) + b / psd^2
  }
  o <- stats::optim(rep(0, ncol(X)), nlp, gr, method = "BFGS", hessian = TRUE)
  S <- solve(o$hessian)
  B <- matrix(stats::rnorm(ndraw * ncol(X)), ndraw) %*% chol(S)
  B <- sweep(B, 2, o$par, "+")
  d1 <- d0 <- data; d1[[trt]] <- 1; d0[[trt]] <- 0
  rhs <- stats::delete.response(stats::terms(formula))
  X1 <- stats::model.matrix(rhs, d1); X0 <- stats::model.matrix(rhs, d0)
  m1 <- colMeans(stats::plogis(X1 %*% t(B)))
  m0 <- colMeans(stats::plogis(X0 %*% t(B)))
  data.frame(m1 = m1, m0 = m0, rd = m1 - m0, log_or = stats::qlogis(m1) - stats::qlogis(m0))
}
