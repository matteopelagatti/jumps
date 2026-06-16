#' Smoothing splines with automatic jump detection via time-warping and
#' Kalman filter-smoother estimation
#'
#' Implements the state-space time-warping model described in the companion
#' paper. The effective inter-observation interval \eqn{\tau_t = \delta_t +
#' \gamma_t} replaces the fixed spacing in both the transition matrix and the
#' process-noise covariance, so jump detection reduces to finding the optimal
#' \eqn{\boldsymbol\gamma} by maximum likelihood.  The analytical gradient of
#' the log-likelihood with respect to \eqn{\gamma_t} is obtained through the
#' Koopman-Shephard (2003) score formula applied to the Kalman smoother output.
#' The optimisation is carried out with the gradient-based CCSAQ algorithm from
#' the NLopt library.
#'
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param lambda positive scalar smoothing parameter (ratio of observation to
#'   process noise variance).
#' @param maxsum maximum sum of \eqn{\gamma_t} (time-warp budget).
#' @param edf logical: if TRUE (default) computes effective degrees of freedom
#'   as the trace of the hat matrix; if FALSE counts non-zero \eqn{\gamma_t}.
#' @param parinit NULL or a numeric vector of length n with starting values for
#'   the optimizer (sigma, gamma_1, ..., gamma_{n-1}); pass the previous
#'   solution for warm-starting.
#' @param last_delta positive scalar: dummy spacing appended after the last
#'   observation (does not affect the fit; default 1).
#' @return list with slots:
#' \itemize{
#'  \item opt: the output of the nloptr optimizer.
#'  \item nobs: number of observations.
#'  \item df: model degrees of freedom.
#'  \item maxsum: the maxsum argument.
#'  \item loglik: Gaussian log-likelihood at the optimum.
#'  \item pars: named vector (sigma, lambda).
#'  \item gamma: estimated time-warp vector of length n-1.
#'  \item ic: vector of information criteria (aic, aicc, bic, hq).
#'  \item smoothed_level: smoothed estimate of f at each x.
#'  \item var_smoothed_level: pointwise variance of the smoothed level.
#'  \item x: the ordered x variable.
#' }
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_kfs <- ssj_kfs(y = faithful$waiting,
#'                   x = faithful$eruptions,
#'                   lambda = 1,
#'                   maxsum = 150)
#' lines(x = sort(faithful$eruptions),
#'       y = of_kfs$smoothed_level, col = "blue", lwd = 2)
#'
#' @export
ssj_kfs <- function(y, x, lambda, maxsum = sd(y, na.rm = TRUE)/mean(diff(x)),
                    edf = TRUE, parinit = NULL, last_delta = 1) {
  n    <- length(y)
  n1   <- n + 1
  nobs <- sum(!is.na(y))
  vy   <- var(y, na.rm = TRUE)
  vdy  <- var(diff(y), na.rm = TRUE)
  sdy  <- sqrt(vy)
  ord  <- order(x)
  x    <- x[ord]
  y    <- y[ord]
  # delta has length n; the last element is the dummy spacing after y[n]
  delta <- c(diff(x), last_delta)

  # Pre-allocate KF/smoother arrays — modified in place by llt_delta
  var_eps      <- numeric(n)
  var_eta      <- numeric(n)
  var_zeta     <- numeric(n)
  cov_eta_zeta <- numeric(n)
  a1  <- rep(y[!is.na(y)][1], n1)
  a2  <- numeric(n1)
  p11 <- rep(vy  * 1e5, n1)
  p12 <- numeric(n1)
  p22 <- rep(vdy * 1e5, n1)
  k1  <- numeric(n)
  k2  <- numeric(n)
  inn <- numeric(n)   # innovations
  fv  <- numeric(n)   # innovation variances
  r1  <- numeric(n1)
  r2  <- numeric(n1)
  n11 <- numeric(n1)
  n12 <- numeric(n1)
  n22 <- numeric(n1)
  e   <- numeric(n1)
  d   <- numeric(n1)
  w   <- numeric(n)
  cnst <- log(2 * pi) / 2

  # Objective (negative log-likelihood / nobs) and its analytical gradient.
  # Parameters: pars[1] = sigma, pars[2:n] = gamma[1:(n-1)].
  obj <- function(pars, wgt = FALSE) {
    sigma <- pars[1]
    gamma <- pars[2:n]                 # length n-1
    sig2  <- sigma^2
    tau   <- delta[1:(n-1)] + gamma    # stretched intervals, length n-1
    tauf  <- c(tau, last_delta)         # with dummy tail, length n

    var_eps[]      <- lambda * sig2
    var_eta[]      <- sig2 * tauf^3 / 3
    var_zeta[]     <- sig2 * tauf
    cov_eta_zeta[] <- sig2 * tauf^2 / 2

    mloglik <- if (wgt) {
      -llt_delta(y, tauf, var_eps, var_eta, var_zeta, cov_eta_zeta,
                 a1, a2, p11, p12, p22, k1, k2, inn, fv,
                 r1, r2, n11, n12, n22, e, d, w)
    } else {
      -llt_delta(y, tauf, var_eps, var_eta, var_zeta, cov_eta_zeta,
                 a1, a2, p11, p12, p22, k1, k2, inn, fv,
                 r1, r2, n11, n12, n22, e, d)
    }

    # Koopman-Shephard (2003) analytical scores via smoother output.
    # r1[-1] = r1[2:(n+1)] stores r_t^(1) for t=1,...,n (shifted by 1).
    rn1  <- r1[-1]^2       - n11[-1]          # r^(1)^2   - N^(11)
    rn2  <- r2[-1]^2       - n22[-1]          # r^(2)^2   - N^(22)
    rn12 <- r1[-1] * r2[-1] - n12[-1]         # r^(1)*r^(2) - N^(12)

    # Score w.r.t. sigma: process terms over t=1..n, observation terms over t=1..n.
    # dQ_t/d(sigma) = 2*sigma * Delta_t  =>  score = sigma * tr((r_t r_t' - N_t) Delta_t)
    process_s  <- sum(rn1 * tauf^3/3 + rn12 * tauf^2 + rn2 * tauf)
    obs_s      <- lambda * sum(e * e - d)   # dH_t/d(sigma) = 2*lambda*sigma
    grad_sigma <- sigma * (process_s + obs_s)

    # Score w.r.t. gamma_t (t=1,...,n-1):
    # dQ_t/d(gamma_t) = sigma^2 * [[tau_t^2, tau_t],[tau_t, 1]]
    grad_gamma <- (sig2 / 2) * (rn1[1:(n-1)] * tau^2 +
                                  2 * rn12[1:(n-1)] * tau +
                                  rn2[1:(n-1)])

    list(
      objective = mloglik / nobs,
      gradient  = c(-grad_sigma, -grad_gamma) / nobs
    )
  }

  # Single inequality constraint: sum(gamma) <= maxsum
  g <- function(pars, wgt) {
    list(
      constraints = sum(pars[2:n]) - maxsum,
      jacobian    = c(0, rep(1, n - 1))
    )
  }

  quasizero <- sdy * 1e-9
  inits <- if (is.null(parinit)) c(sdy / 10, rep(0, n - 1)) else parinit
  lb    <- c(quasizero, rep(0, n - 1))
  inits <- pmax(inits, lb)

  opt <- nloptr::nloptr(
    x0          = inits,
    eval_f      = obj,
    lb          = lb,
    eval_g_ineq = g,
    opts        = list(algorithm         = "NLOPT_LD_CCSAQ",
                       xtol_rel          = 1e-5,
                       check_derivatives = FALSE,
                       maxeval           = 2000),
    wgt = FALSE
  )

  if (edf) {
    obj(opt$solution, TRUE)   # populates w via llt_delta side-effect
    df <- sum(w)
  } else {
    df <- 1 + sum(opt$solution[2:n] > quasizero)
  }
  loglik <- -nobs * (opt$objective + cnst)

  list(
    opt                = opt,
    nobs               = n,
    df                 = df,
    maxsum             = maxsum,
    loglik             = loglik,
    pars               = c(sigma = opt$solution[1], lambda = lambda),
    gamma              = opt$solution[2:n],
    ic                 = c(aic  =  2 * (df - loglik),
                           aicc =  2 * (df * n / (n - df - 1) - loglik),
                           bic  =  df * log(n) - 2 * loglik,
                           hq   =  2 * (log(log(n)) * df - loglik)),
    smoothed_level     = (a1 + p11 * r1 + p12 * r2)[-(n + 1)],
    var_smoothed_level = (p11 - p11^2 * n11 - 2 * p11 * p12 * n12 -
                            p12^2 * n22)[-(n + 1)],
    x                  = x
  )
}


#' Automatic selection of the optimal time-warping smoothing spline with jumps
#' using the Kalman filter-smoother
#'
#' Runs \code{ssj_kfs} over a grid of \code{maxsum} values and returns the
#' solution with the lowest selected information criterion.  Each grid point is
#' warm-started from the previous optimum.  The search stops early when the IC
#' increases for three consecutive grid values.
#'
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param lambda smoothing constant (positive scalar).
#' @param grid numeric vector of maxsum values to search; default spans 0 to
#'   10*sd(y)/mean(diff(x)) in steps of sd(y)/(10*mean(diff(x))).
#' @param ic string: information criterion for selection ("bic", "hq",
#'   "aic", "aicc").
#' @param edf logical: passed to \code{ssj_kfs} (default TRUE).
#' @param last_delta numeric scalar: passed to \code{ssj_kfs} (default 1).
#'
#' @returns The output of \code{ssj_kfs} corresponding to the best IC.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_kfs <- auto_ssj_kfs(y = faithful$waiting,
#'                        x = faithful$eruptions,
#'                        lambda = 1)
#' lines(x = sort(faithful$eruptions),
#'       y = of_kfs$smoothed_level, col = "blue", lwd = 2)
#'
#' @export
auto_ssj_kfs <- function(y, x, lambda,
                         grid = seq(0,
                                    sd(y, na.rm = TRUE) / mean(diff(x)) * 10,
                                    sd(y, na.rm = TRUE) / mean(diff(x)) / 10),
                         ic = c("bic", "hq", "aic", "aicc"),
                         edf = TRUE, last_delta = 1) {
  ic          <- match.arg(ic)
  last_ic     <- Inf
  n_increases <- 0L
  parinit     <- NULL
  for (M in grid) {
    out <- ssj_kfs(y = y, x = x, lambda = lambda,
                   maxsum = M, edf = edf, parinit = parinit,
                   last_delta = last_delta)
    parinit    <- c(out$pars["sigma"], out$gamma)   # warm start
    current_ic <- switch(ic,
                         bic  = out$ic["bic"],
                         hq   = out$ic["hq"],
                         aic  = out$ic["aic"],
                         aicc = out$ic["aicc"])
    if (current_ic < last_ic) {
      best        <- out
      last_ic     <- current_ic
      n_increases <- 0L
    } else {
      n_increases <- n_increases + 1L
      if (n_increases >= 3L) break
    }
  }
  best
}
