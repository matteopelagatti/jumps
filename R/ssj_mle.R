#' Smoothing splines with automatic jump detection via time-warping, with
#' the smoothing parameter lambda estimated by maximum likelihood
#'
#' Identical to \code{ssj_kfs} except that \code{lambda} is no longer a fixed
#' input: it is a free parameter optimized jointly with sigma and gamma. Since
#' lambda enters only the observation variance \eqn{H_t = \lambda\sigma^2} and
#' never the transition matrix, its score does not need the extra correction
#' that gamma requires -- the plain Koopman-Shephard formula applies directly:
#' \eqn{\partial\ell/\partial\lambda = (\sigma^2/2)\sum_t(e_t^2-d_t)}. This
#' avoids having to grid-search lambda (which would suffer from the same
#' degrees-of-freedom problem as maxsum) on top of the maxsum grid search.
#' Internally lambda is optimized via \eqn{\mu=\log\lambda}, which guarantees
#' lambda > 0 without a box constraint and avoids the scale instability seen
#' when optimizing the raw lambda directly (its sensible range can span
#' orders of magnitude, like a variance ratio); the score is a one-line
#' chain-rule correction, \eqn{\partial\ell/\partial\mu =
#' \lambda\,\partial\ell/\partial\lambda}. \code{parinit} and the returned
#' \code{pars} are unaffected and stay in terms of the raw lambda. When
#' \code{parinit} is not supplied, the optimizer's success turns out to
#' depend on the path it traverses, not just on how close the starting
#' lambda looks to the eventual optimum -- so rather than rely on a single
#' formula, \code{ssj_mle} tries a handful of starting lambdas spanning
#' several orders of magnitude (anchored at a robust, jump-resistant noise
#' estimate via \code{mad(diff(y))}) and keeps whichever converges to the
#' highest likelihood.
#'
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param maxsum maximum sum of \eqn{\gamma_t} (time-warp budget).
#' @param edf logical: if TRUE (default) computes effective degrees of freedom
#'   as the trace of the hat matrix; if FALSE counts non-zero \eqn{\gamma_t}.
#' @param parinit NULL or a numeric vector of length n+1 with starting values
#'   for the optimizer (sigma, lambda, gamma_1, ..., gamma_{n-1}); pass the
#'   previous solution for warm-starting.
#' @param last_delta positive scalar: dummy spacing appended after the last
#'   observation (does not affect the fit; default 1).
#' @param ebic_xi numeric scalar in [0,1]: strength of the Extended-BIC penalty
#'   for searching over the n-1 candidate jump locations (default 1, the
#'   Chen and Chen 2008 recommendation for sparse selection; 0 reduces ebic to
#'   plain bic).
#' @return list with slots:
#' \itemize{
#'  \item opt: the output of the nloptr optimizer.
#'  \item nobs: number of observations.
#'  \item df: model degrees of freedom.
#'  \item maxsum: the maxsum argument.
#'  \item loglik: Gaussian log-likelihood at the optimum.
#'  \item pars: named vector (sigma, lambda); lambda is now estimated.
#'  \item gamma: estimated time-warp vector of length n-1.
#'  \item ic: vector of information criteria (aic, aicc, bic, hq, ebic).
#'  \item smoothed_level: smoothed estimate of f at each x.
#'  \item var_smoothed_level: pointwise variance of the smoothed level.
#'  \item x: the ordered x variable.
#' }
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_mle <- ssj_mle(y = faithful$waiting,
#'                   x = faithful$eruptions,
#'                   maxsum = 150)
#' lines(x = sort(faithful$eruptions),
#'       y = of_mle$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
ssj_mle <- function(y, x, maxsum = sd(y, na.rm = TRUE)/mean(diff(x)),
                    edf = TRUE, parinit = NULL, last_delta = 1, ebic_xi = 1) {
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

  # Normalise to make the problem scale-independent: x -> [0,1], y -> unit
  # variance.  Units: sigma ~ [y]/[x]^{3/2}, so sigma_norm = sigma_orig *
  # x_range^{3/2} / y_scale; lambda ~ [x]^3, so lambda_norm = lambda_orig /
  # x_range^3; gamma_norm = gamma_orig / x_range.  All back-transformations are
  # applied to the returned quantities; the optimiser works in normalised units.
  # Consequence: lambda_init_candidates can be a fixed grid independent of the
  # data scale, eliminating the instability of the 1/dbar^3 formula.
  x_range    <- x[n] - x[1L]   # x is already sorted
  x_min_s    <- x[1L]
  y_scale    <- sdy
  if (x_range < .Machine$double.eps) stop("'x' has zero range")
  x          <- (x - x_min_s) / x_range
  y          <- y / y_scale
  delta      <- delta / x_range
  last_delta <- last_delta / x_range  # normalise the dummy tail spacing
  maxsum_orig <- maxsum
  maxsum      <- maxsum / x_range     # normalise the budget constraint
  # Update statistics to normalised coordinates
  vy  <- 1.0           # var(y/y_scale) = 1 by construction
  sdy <- 1.0
  vdy <- vdy / y_scale^2

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
  # Parameters: pars[1] = sigma, pars[2] = lambda, pars[3:(n+1)] = gamma[1:(n-1)].
  obj <- function(pars, wgt = FALSE) {
    sigma  <- pars[1]
    lambda <- exp(pars[2])              # pars[2] = mu = log(lambda): guarantees
                                          # lambda > 0 with no box constraint, and
                                          # avoids the scale blow-up seen when
                                          # optimizing the raw (additive) lambda.
    gamma  <- pars[3:(n + 1)]           # length n-1
    sig2   <- sigma^2
    tau    <- delta[1:(n-1)] + gamma    # stretched intervals, length n-1
    tauf   <- c(tau, last_delta)         # with dummy tail, length n

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
    rn1  <- r1[-1]^2       - n11[-1]
    rn2  <- r2[-1]^2       - n22[-1]
    rn12 <- r1[-1] * r2[-1] - n12[-1]

    # H_t = lambda*sigma^2 enters ONLY the observation variance, never the
    # transition matrix, so both sigma's and lambda's scores follow the plain
    # Koopman-Shephard formula with no extra correction term -- unlike gamma.
    process_s   <- sum(rn1 * tauf^3/3 + rn12 * tauf^2 + rn2 * tauf)
    obs_s       <- sum(e * e - d)                 # raw, no lambda factor
    grad_sigma  <- sigma * (process_s + lambda * obs_s)   # dH_t/d(sigma)  = 2*lambda*sigma
    grad_lambda <- (sig2 / 2) * obs_s                      # dH_t/d(lambda) = sigma^2
    grad_mu     <- lambda * grad_lambda   # chain rule: d(loglik)/d(mu) = lambda * d(loglik)/d(lambda)
                                            # (verified numerically against numDeriv, incl. at lambda ~ 2e4)

    # Score w.r.t. gamma_t: Qterm is the plain Koopman-Shephard piece (Q_t
    # dependence); Tterm is the extra piece from gamma_t also entering the
    # transition matrix T_t -- see ssj_kfs.R for the full derivation.
    Qterm <- (sig2 / 2) * (rn1[1:(n-1)] * tau^2 +
                             2 * rn12[1:(n-1)] * tau +
                             rn2[1:(n-1)])

    kk        <- 1 - k1[1:(n-1)]
    alphahat2 <- a2[1:(n-1)] + p12[1:(n-1)] * r1[1:(n-1)] + p22[1:(n-1)] * r2[1:(n-1)]
    PL21      <- p12[1:(n-1)] * (kk * n11[2:n] - k2[1:(n-1)] * n12[2:n]) +
                 p22[1:(n-1)] * (tau * n11[2:n] + n12[2:n])
    Tterm     <- alphahat2 * r1[2:n] - PL21

    grad_gamma <- Qterm + Tterm

    list(
      objective = mloglik / nobs,
      gradient  = c(-grad_sigma, -grad_mu, -grad_gamma) / nobs
    )
  }

  # Single inequality constraint: sum(gamma) <= maxsum
  g <- function(pars, wgt) {
    list(
      constraints = sum(pars[3:(n + 1)]) - maxsum,
      jacobian    = c(0, 0, rep(1, n - 1))
    )
  }

  quasizero <- sdy * 1e-9
  lb        <- c(quasizero, -Inf, rep(0, n - 1))

  # No single closed-form lambda guess proved robust across datasets: the
  # optimizer's success turns out to depend on the *path* it traverses, not
  # just on how close the starting lambda looks to the eventual optimum --
  # a deliberately "wrong" (very large) starting lambda sometimes converges
  # correctly precisely because traversing down from it passes through
  # regions where gamma's benefit becomes visible, while starting already
  # close to the optimal lambda can let the optimizer settle near gamma = 0
  # without ever exploring further. Rather than chase a single formula, try
  # a handful of starting lambdas spanning several orders of magnitude and
  # keep whichever converges to the best (highest) likelihood -- the same
  # "best of several attempts" principle already used below for the
  # algorithm choice. This only applies when no explicit parinit is given;
  # a user- or auto_ssj_mle-supplied parinit is used as-is, with no
  # multi-start, since it already encodes a specific, informed choice (e.g.
  # warm-started from the previous grid point).
  # After normalising x -> [0,1] and y -> unit variance, lambda_norm lives in
  # a dataset-independent range.  A fixed six-point grid spanning 1e-4 to 1e2
  # covers the practically relevant region without a fragile 1/dbar^3 formula.
  lambda_init_candidates <- function() c(1e-4, 1e-2, 0.1, 1, 10, 100)

  run_nloptr <- function(opts, x0) {
    nloptr::nloptr(x0 = x0, eval_f = obj, lb = lb, eval_g_ineq = g,
                   opts = opts, wgt = FALSE)
  }

  # AUGLAG wrapping LBFGS converges far faster than CCSAQ on this problem when
  # it works, but it has its own failure mode: it can report "converged"
  # (ftol_reached) sitting at gamma = 0 with a large unused budget and a
  # strongly improving direction still available -- a status-code check alone
  # would not catch that. So in addition to the status code, check the KKT
  # condition directly: with the budget constraint slack, every gamma_t stuck
  # at its lower bound (0) must have non-positive score, or the claimed
  # optimum is spurious and we fall back to CCSAQ.
  looks_unconverged <- function(sol) {
    gamma       <- sol[3:(n + 1)]
    budget_left <- maxsum - sum(gamma)
    if (budget_left < 1e-8 * max(maxsum, 1)) return(FALSE)  # budget fully used
    inactive <- gamma <= 1e-8 * max(maxsum, 1)
    if (!any(inactive)) return(FALSE)
    grad_gamma <- -obj(sol)$gradient[3:(n + 1)] * nobs   # raw score d(loglik)/d(gamma)
    scale <- max(abs(grad_gamma), 1e-8)
    max(grad_gamma[inactive]) > 1e-6 * scale
  }

  solve_with_fallback <- function(x0) {
    opt <- run_nloptr(list(algorithm         = "NLOPT_LD_AUGLAG",
                           xtol_rel          = 1e-5,
                           check_derivatives = FALSE,
                           maxeval           = 2000,
                           local_opts        = list(algorithm = "NLOPT_LD_LBFGS",
                                                    xtol_rel  = 1e-6)),
                      x0)
    if (opt$status < 0 || opt$status == 5 || looks_unconverged(opt$solution)) {
      opt <- run_nloptr(list(algorithm         = "NLOPT_LD_CCSAQ",
                             xtol_rel          = 1e-5,
                             check_derivatives = FALSE,
                             maxeval           = 2000),
                        x0)
    }
    opt
  }

  make_x0 <- function(lambda0) {
    x0 <- c(sdy / 10, lambda0, rep(0, n - 1))
    x0[2] <- log(max(x0[2], quasizero))
    pmax(x0, lb)
  }

  # parinit is accepted in ORIGINAL units and converted here to normalised
  # units so callers (including auto_ssj_mle) can always work in the natural
  # scale of the data.  The returned pars / gamma are back-transformed below.
  if (is.null(parinit)) {
    attempts <- lapply(lambda_init_candidates(), function(lam0) solve_with_fallback(make_x0(lam0)))
    opt <- attempts[[which.min(vapply(attempts, function(a) a$objective, numeric(1)))]]
  } else {
    x0      <- parinit
    x0[1L]             <- x0[1L] * x_range^1.5 / y_scale   # sigma:  orig -> norm
    x0[2L]             <- x0[2L] / x_range^3               # lambda: orig -> norm (pre-log)
    x0[3L:(n + 1L)]   <- x0[3L:(n + 1L)] / x_range        # gamma:  orig -> norm
    x0[2L] <- log(max(x0[2L], quasizero))
    x0 <- pmax(x0, lb)
    opt <- solve_with_fallback(x0)
  }

  if (edf) {
    obj(opt$solution, TRUE)   # populates w, r1, r2, n11, n12, n22 via llt_delta side-effect
    df <- sum(w)
    # Smoothed level and slope disturbances: E[eta_t | Y] = Q_t * [r1_{t+1}; r2_{t+1}].
    # var_eta, var_zeta, cov_eta_zeta are already set by the llt_delta call above.
    # Length n (index t = 1,...,n); only t = 1,...,n-1 are meaningful for gamma.
    eta_hat  <- var_eta * r1[-1] + cov_eta_zeta * r2[-1]
    zeta_hat <- cov_eta_zeta * r1[-1] + var_zeta * r2[-1]
    # Variance of smoothed disturbances: Var(eta_hat_t) = Q_t - Q_t N_{t+1} Q_t
    var_eta_hat  <- pmax(var_eta  - (var_eta^2        * n11[-1] +
                                     2*var_eta*cov_eta_zeta * n12[-1] +
                                     cov_eta_zeta^2   * n22[-1]), 0)
    var_zeta_hat <- pmax(var_zeta - (cov_eta_zeta^2   * n11[-1] +
                                     2*cov_eta_zeta*var_zeta * n12[-1] +
                                     var_zeta^2       * n22[-1]), 0)
    eps <- .Machine$double.eps
    std_eta_hat  <- eta_hat  / sqrt(var_eta_hat  + eps)
    std_zeta_hat <- zeta_hat / sqrt(var_zeta_hat + eps)
  } else {
    df <- 2 + sum(opt$solution[3:(n + 1)] > quasizero)
    eta_hat <- zeta_hat <- std_eta_hat <- std_zeta_hat <- NULL
  }
  loglik <- -nobs * (opt$objective + cnst)

  # Extended BIC (Chen and Chen, 2008): bic plus 2*xi*k*log(p), where k is the
  # number of active jump locations and p = n-1 the number of candidates.
  k_active <- sum(opt$solution[3:(n + 1)] > quasizero)
  bic_val  <- df * log(n) - 2 * loglik
  ebic_val <- bic_val + 2 * ebic_xi * k_active * log(n - 1)

  # Back-transform all returned quantities to original units.
  list(
    opt                = opt,
    nobs               = n,
    df                 = df,
    maxsum             = maxsum_orig,
    loglik             = loglik,
    pars               = c(sigma  = opt$solution[1L] * y_scale / x_range^1.5,
                           lambda = exp(opt$solution[2L]) * x_range^3),
    gamma              = opt$solution[3L:(n + 1L)] * x_range,
    ic                 = c(aic  =  2 * (df - loglik),
                           aicc =  2 * (df * n / (n - df - 1) - loglik),
                           bic  =  bic_val,
                           hq   =  2 * (log(log(n)) * df - loglik),
                           ebic =  ebic_val),
    smoothed_level     = (a1 + p11 * r1 + p12 * r2)[-(n + 1L)] * y_scale,
    var_smoothed_level = (p11 - p11^2 * n11 - 2 * p11 * p12 * n12 -
                            p12^2 * n22)[-(n + 1L)] * y_scale^2,
    dist_level         = eta_hat  * y_scale,          # level disturbance in orig y units
    dist_slope         = zeta_hat * y_scale / x_range, # slope disturbance in orig y/x units
    std_dist_level     = std_eta_hat,                  # dimensionless, no rescaling
    std_dist_slope     = std_zeta_hat,
    x                  = x * x_range + x_min_s        # back to original x
  )
}


#' Automatic selection of the optimal time-warping smoothing spline with jumps
#' using maximum-likelihood estimation of sigma, lambda and the time-warps
#'
#' Runs \code{ssj_mle} over a grid of \code{maxsum} values and returns the
#' solution with the lowest selected information criterion. Unlike
#' \code{auto_ssj_kfs}/\code{auto_ssj_atw}, lambda is not part of the grid: it
#' is estimated by maximum likelihood inside every \code{ssj_mle} call, so
#' only maxsum needs to be searched. Each grid point is warm-started from the
#' previous optimum (sigma, lambda and gamma). The search stops early when the
#' IC increases for three consecutive grid values.
#'
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param grid numeric vector of maxsum values to search; default spans 0 to
#'   10*sd(y)/mean(diff(x)) in steps of sd(y)/(10*mean(diff(x))).
#' @param ic string: information criterion for selection. Default "ebic"
#'   (recommended when jumps are expected to be rare relative to the n-1
#'   candidate locations -- see \code{ebic_xi}); "bic", "hq", "aic" and
#'   "aicc" are also available.
#' @param edf logical: passed to \code{ssj_mle} (default TRUE).
#' @param last_delta numeric scalar: passed to \code{ssj_mle} (default 1).
#' @param ebic_xi numeric scalar in [0,1] passed to \code{ssj_mle}: strength
#'   of the Extended-BIC search penalty (default 1; only relevant when
#'   \code{ic = "ebic"}).
#'
#' @returns The output of \code{ssj_mle} corresponding to the best IC.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_mle <- auto_ssj_mle(y = faithful$waiting,
#'                        x = faithful$eruptions)
#' lines(x = sort(faithful$eruptions),
#'       y = of_mle$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
auto_ssj_mle <- function(y, x,
                         grid = seq(0, max(x) - min(x), length.out = 21L),
                         ic = c("ebic", "bic", "hq", "aic", "aicc"),
                         edf = TRUE, last_delta = 1, ebic_xi = 1) {
  ic <- match.arg(ic)
  n  <- length(y)

  # Stage 1: fit the pure smooth (maxsum = 0).  With no gamma degrees of
  # freedom the optimisation landscape is unimodal and sigma/lambda are
  # identified cleanly.  The smoothed level disturbances E[eta_t | Y] are
  # large at locations where the smooth model is "surprised" by a sudden
  # level change — i.e. at genuine discontinuities.
  fit0   <- ssj_mle(y = y, x = x, maxsum = 0, edf = TRUE,
                    last_delta = last_delta, ebic_xi = ebic_xi)
  sigma0 <- fit0$pars["sigma"]
  lam0   <- fit0$pars["lambda"]
  # abs level disturbances for t = 1, ..., n-1  (t = n is the dummy tail)
  # dist_level is returned in original y units; proportional allocation below
  # gives gamma init values that sum to M (in original x units), which ssj_mle
  # normalises internally.
  eta_abs <- abs(fit0$dist_level[seq_len(n - 1L)])

  # Allocate the maxsum budget proportionally to the disturbance magnitudes,
  # giving the optimizer a warm start at the right *locations*.
  init_gamma <- function(M) {
    s <- sum(eta_abs)
    if (s < .Machine$double.eps) return(rep(0, n - 1L))
    eta_abs / s * M
  }

  last_ic     <- Inf
  n_increases <- 0L
  best        <- fit0

  for (M in grid) {
    if (M == 0) {
      out <- fit0
    } else {
      parinit <- c(sigma0, lam0, init_gamma(M))
      out <- ssj_mle(y = y, x = x,
                     maxsum = M, edf = edf, parinit = parinit,
                     last_delta = last_delta, ebic_xi = ebic_xi)
    }
    current_ic <- switch(ic,
                         bic  = out$ic["bic"],
                         ebic = out$ic["ebic"],
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
