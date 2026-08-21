#' Smoothing splines with automatic jump detection via time warping and
#' Kalman filter-smoother estimation
#'
#' Implements the state-space time-warping model described in the companion
#' paper. The effective inter-observation interval \eqn{\tau_t = \delta_t +
#' \gamma_t} replaces the fixed spacing in both the transition matrix and the
#' process-noise covariance, so that jump detection reduces to finding the
#' optimal \eqn{\gamma} by maximum likelihood.
#'
#' The analytical gradient of the log-likelihood with respect to \eqn{\gamma_t}
#' is obtained through the Koopman-Shephard (2003) score formula applied to the
#' Kalman smoother output, plus an extra term accounting for \eqn{\gamma_t}
#' entering the transition matrix. The optimisation is carried out with the
#' gradient-based AUGLAG/LBFGS algorithms of the NLopt library, with a fallback
#' to CCSAQ whenever the KKT conditions are not satisfied at the reported
#' optimum.
#'
#'
#' Tied values of \code{x} are collapsed to their weighted mean before fitting,
#' and the fitted values are returned at the vector of unique, sorted \code{x}
#' available in the \code{x} slot of the output. This keeps the problem well
#' posed: with ties the interval between two replicates has zero length, and a
#' warp placed there would be a discontinuity between two observations recorded
#' at the same \code{x}, which is not identified. The replicate counts become
#' observation weights automatically, so no information is lost.
#' @param y numeric vector with the response (missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param lambda positive scalar smoothing parameter (ratio of observation to
#'   process noise variance).
#' @param maxsum maximum sum of the \eqn{\gamma_t}, that is the warping budget
#'   \eqn{M} of the constraint \eqn{\sum_t \gamma_t \le M}. It is expressed in
#'   the units of \code{x} and defaults to the range of \code{x}, the same
#'   upper end as the grid searched by \code{\link{auto_ssj_kfs}}.
#' @param weights optional numeric vector of positive observation weights, one
#'   per element of \code{y}. The observation variance becomes
#'   \eqn{H_t = \lambda\sigma^2/w_t}, so a point carrying weight \eqn{w_t} is
#'   treated as the mean of \eqn{w_t} independent observations. This is what
#'   makes it correct to collapse tied values of \code{x} to their mean and
#'   pass the number of replicates as the weight: the Gaussian likelihood
#'   factorises, so the weighted fit on the collapsed data is identical to the
#'   fit on the raw replicates. The default, \code{NULL}, means unit weights.
#' @param edf logical: if \code{TRUE} (default) the effective degrees of freedom
#'   are computed as the trace of the hat matrix; if \code{FALSE} the non-zero
#'   \eqn{\gamma_t} are counted.
#' @param parinit \code{NULL} or a numeric vector of length \code{n} with
#'   starting values for the optimizer, namely \eqn{(\sigma, \gamma_1, \ldots,
#'   \gamma_{n-1})}; pass the previous solution for warm starting.
#' @param last_delta positive scalar: dummy spacing appended after the last
#'   observation (it does not affect the fit; default 1).
#' @param ebic_xi numeric scalar in \code{[0,1]}: strength of the extended-BIC
#'   penalty for searching over the \code{n-1} candidate jump locations
#'   (default 0.15; 0 reduces \code{ebic} to plain \code{bic}). The default is
#'   calibrated rather than taken from theory; see the note below.
#'
#' @section Choosing \code{ebic_xi}:
#' The default is the value that minimises worst-case regret in root mean
#' squared error over the 48-cell design of
#' \code{inst/scripts/xi_calibration.R}. It is lower than the 0.3 used by
#' \code{\link{ssj_mle}}, and the reason is not that the models differ -- at a
#' given \code{lambda} and budget they are the same fit -- but that the
#' procedures do. \code{\link{auto_ssj_kfs}} holds \code{lambda} fixed across
#' the budget grid, whereas \code{\link{auto_ssj_mle}} re-estimates it at every
#' point; being able to re-tune the smoothness as the budget grows gives more
#' scope to overfit, and so calls for a heavier penalty.
#'
#' The value sits below the Chen and Chen (2008) consistency floor of about 0.5
#' because the plain BIC term already charges for each activated warp through
#' \code{df}, at about \eqn{0.34\log n} per declared jump for this estimator --
#' an implicit charge of roughly 0.17 on the \eqn{\xi} scale.
#'
#' @returns An object of class \code{ssj} (a list, so that all the slots below
#'   remain directly accessible) with a \code{\link[=print.ssj]{print}} and a
#'   \code{\link[=plot.ssj]{plot}} method. Its slots are:
#' \itemize{
#'  \item \code{opt}: the output of the \code{nloptr} optimizer.
#'  \item \code{nobs}: number of observations.
#'  \item \code{df}: model's degrees of freedom.
#'  \item \code{maxsum}: the \code{maxsum} argument.
#'  \item \code{loglik}: Gaussian log-likelihood at the optimum, on the scale
#'    of \code{y} as supplied. Note that the diffuse initial conditions are
#'    approximated by a large proper prior, so the absolute value carries an
#'    additive constant set by that convention; differences across
#'    \code{maxsum} for one data set, which is what the information criteria
#'    compare, are unaffected.
#'  \item \code{pars}: named vector (sigma, lambda).
#'  \item \code{gamma}: estimated time-warp vector of length \code{n-1}.
#'  \item \code{ic}: vector of information criteria (aic, aicc, bic, hq, ebic).
#'  \item \code{smoothed_level}: smoothed estimate of \eqn{f} at each \code{x}.
#'  \item \code{var_smoothed_level}: pointwise variance of the smoothed level.
#'  \item \code{x}: the ordered \code{x} variable.
#' }
#'
#' @references
#' Chen, J. and Chen, Z. (2008). Extended Bayesian information criteria for
#' model selection with large model spaces. \emph{Biometrika}, \strong{95}(3),
#' 759--771. \doi{10.1093/biomet/asn034}
#'
#' Koopman, S. J. and Shephard, N. (1992). Exact score for time series models in
#' state space form. \emph{Biometrika}, \strong{79}(4), 823--826.
#' \doi{10.1093/biomet/79.4.823}
#'
#' @seealso \code{\link{auto_ssj_kfs}}, \code{\link{ssj_mle}},
#'   \code{\link{ssj_atw}}, \code{\link{ssj_fix}}.
#'
#' @examples
#' # 'maxsum' is a warping budget in the units of x; the range of eruptions
#' # is about 3.5, so a budget of 1 is already generous.
#' plot(faithful$eruptions, faithful$waiting)
#' of_kfs <- ssj_kfs(y = faithful$waiting,
#'                   x = faithful$eruptions,
#'                   lambda = 1,
#'                   maxsum = 1)
#' lines(x = of_kfs$x, y = of_kfs$smoothed_level, col = "blue", lwd = 2)
#'
#' @export
ssj_kfs <- function(y, x, lambda, maxsum = max(x) - min(x), weights = NULL,
                    edf = TRUE, parinit = NULL, last_delta = 1, ebic_xi = 0.15) {
  wt <- ssj_check_weights(weights, length(y))
  y  <- as.numeric(y); x <- as.numeric(x)
  ord  <- order(x)
  x    <- x[ord]
  y    <- y[ord]
  wt   <- wt[ord]
  # Tied values of x are collapsed to their weighted mean.  Doing it here rather
  # than leaving it to the caller keeps the problem well posed: with ties the
  # intervals between replicates have zero length, and a warp placed on one of
  # them is a discontinuity between two observations recorded at the same x,
  # which is not identified.
  cl_  <- ssj_collapse_ties(y, x, wt)
  y    <- cl_$y; x <- cl_$x; wt <- cl_$wt
  n    <- length(y)
  n1   <- n + 1
  nobs <- sum(!is.na(y))
  vy   <- var(y, na.rm = TRUE)
  sdy  <- sqrt(vy)
  wf   <- c(wt, 1)   # padded: e and d carry a trailing unused element
  # NB: computed AFTER sorting.  var(diff(y)) is not invariant to the order in
  # which the observations are supplied, and it scales the diffuse prior on the
  # slope (p22 below), so taking it on the unsorted y would make the fit depend
  # on the row order of the input.
  vdy  <- var(diff(y), na.rm = TRUE)
  # Force the default of maxsum NOW, while x is still on its original scale.
  # The default refers to x, and R evaluates it lazily at first use; without
  # this, the promise would be forced below, after x has been rescaled to
  # [0, 1], and the budget would silently become 1 instead of the range of x.
  force(maxsum)
  # delta has length n; the last element is the dummy spacing after y[n]
  delta <- c(diff(x), last_delta)

  # Normalise: x -> [0,1], y -> unit variance (same logic as ssj_mle).
  # lambda ~ [x]^3 => lambda_norm = lambda / x_range^3.
  x_range    <- x[n] - x[1L]
  x_min_s    <- x[1L]
  y_scale    <- sdy
  if (x_range < .Machine$double.eps) stop("'x' has zero range")
  x          <- (x - x_min_s) / x_range
  y          <- y / y_scale
  delta      <- delta / x_range
  last_delta <- last_delta / x_range
  lambda_orig <- lambda
  lambda      <- lambda / x_range^3
  maxsum_orig <- maxsum
  maxsum      <- maxsum / x_range
  vy  <- 1.0; sdy <- 1.0; vdy <- vdy / y_scale^2

  # Normalise parinit = [sigma, gamma_1,...,gamma_{n-1}] to scaled units
  ssj_check_parinit(parinit, n, "(sigma, gamma_1, ..., gamma_{n-1})",
                    n, cl_$collapsed)
  if (!is.null(parinit)) {
    parinit[1L]   <- parinit[1L] * x_range^1.5 / y_scale   # sigma:  orig -> norm
    parinit[2L:n] <- parinit[2L:n] / x_range                # gamma:  orig -> norm
  }

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

    var_eps[]      <- lambda * sig2 / wt
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
    obs_s      <- lambda * sum((e * e - d) / wf)  # dH_t/d(sigma) = 2*lambda*sigma/w
    grad_sigma <- sigma * (process_s + obs_s)

    # Score w.r.t. gamma_t (t=1,...,n-1). Unlike sigma, gamma_t enters BOTH the
    # process noise Q_t = sigma^2*Delta_t(tau_t) AND the transition matrix
    # T_t = [[1,tau_t],[0,1]]; the plain Koopman-Shephard formula below (Qterm)
    # only accounts for the Q_t-dependence. The T_t-dependence needs the extra
    # smoothed cross-moment Cov(alpha_t,eta_t|y) = -P_t L_t' N_t Q_t, derived via
    # the Fisher/Louis identity and verified against an exact (brute-force
    # Gaussian) smoother; the Q_t / Q_t^{-1} factors cancel by the cyclic trace
    # property, leaving Tterm purely in terms of quantities llt_delta already
    # computes (no new recursion needed).
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
  lb        <- c(quasizero, rep(0, n - 1))

  # Starting value for sigma.  Because lambda is fixed here, sigma is a pure
  # scale: multiplying sigma^2 by c multiplies H_t and Q_t, hence every
  # innovation variance F_t, by c.  The profile log-likelihood in c is then
  # -(1/2) sum(log(c F_t) + v_t^2/(c F_t)), maximised at c = mean(v_t^2/F_t)
  # evaluated once at sigma = 1.  A single filter pass therefore lands
  # essentially on the optimum.
  #
  # This matters.  The previous fixed start (sdy/10, i.e. 0.1 in the units
  # used here) could sit orders of magnitude away: with a small lambda the
  # optimal sigma is large -- 37 on the test case below -- and the optimizer
  # then diverged, returning sigma ~ 1e6 and a log-likelihood 800 units worse
  # than ssj_mle() reached on the same data at the same lambda.  It did so in
  # 37 of 40 test replications, and auto_ssj_kfs() propagated the failure
  # along the whole grid through its warm start.
  profile_sigma <- function() {
    var_eps[]      <- lambda / wt          # sigma = 1
    var_eta[]      <- delta^3 / 3
    var_zeta[]     <- delta
    cov_eta_zeta[] <- delta^2 / 2
    llt_delta(y, delta, var_eps, var_eta, var_zeta, cov_eta_zeta,
              a1, a2, p11, p12, p22, k1, k2, inn, fv,
              r1, r2, n11, n12, n22, e, d)
    ok <- is.finite(fv) & fv > 0 & is.finite(inn)
    if (!any(ok)) return(sdy / 10)
    max(sqrt(mean(inn[ok]^2 / fv[ok])), quasizero)
  }

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
    gamma       <- sol[2:n]
    budget_left <- maxsum - sum(gamma)
    if (budget_left < 1e-8 * max(maxsum, 1)) return(FALSE)  # budget fully used
    inactive <- gamma <= 1e-8 * max(maxsum, 1)
    if (!any(inactive)) return(FALSE)
    grad_gamma <- -obj(sol)$gradient[2:n] * nobs   # raw score d(loglik)/d(gamma)
    scale <- max(abs(grad_gamma), 1e-8)
    max(grad_gamma[inactive]) > 1e-6 * scale
  }

  auglag <- list(algorithm         = "NLOPT_LD_AUGLAG",
                 xtol_rel          = 1e-5,
                 check_derivatives = FALSE,
                 maxeval           = 2000,
                 local_opts        = list(algorithm = "NLOPT_LD_LBFGS",
                                          xtol_rel  = 1e-6))
  ccsaq <- list(algorithm         = "NLOPT_LD_CCSAQ",
                xtol_rel          = 1e-5,
                check_derivatives = FALSE,
                maxeval           = 2000)

  if (is.null(parinit)) {
    # Multi-start around the profile value, in the same spirit as the grid of
    # lambda starts in ssj_mle(): the profile solves the sigma direction, but
    # with maxsum > 0 the gamma surface is not convex and the starting scale
    # still matters.  Keep whichever start reaches the highest likelihood.
    s0       <- profile_sigma()
    starts   <- lapply(s0 * c(1, 0.1, 10),
                       function(s) pmax(c(s, rep(0, n - 1)), lb))
    attempts <- lapply(starts, function(x0) run_nloptr(auglag, x0))
    opt      <- attempts[[which.min(vapply(attempts,
                  function(a) if (is.finite(a$objective)) a$objective else Inf,
                  numeric(1)))]]
    inits    <- starts[[1L]]
  } else {
    inits <- pmax(parinit, lb)
    opt   <- run_nloptr(auglag, inits)
  }

  if (opt$status < 0 || opt$status == 5 || looks_unconverged(opt$solution)) {
    alt <- run_nloptr(ccsaq, inits)
    if (is.finite(alt$objective) && alt$objective < opt$objective) opt <- alt
  }

  if (edf) {
    obj(opt$solution, TRUE)   # populates w via llt_delta side-effect
    df <- sum(w)
  } else {
    df <- 1 + sum(opt$solution[2:n] > quasizero)
  }
  # Back to the scale of the data as supplied: y was divided by y_scale on
  # entry, so the density of the rescaled copy differs from that of the
  # original by one factor of y_scale per observation entering the likelihood.
  # Unlike ssj_mle(), this fitter does not restore the within-group scatter, so
  # the likelihood is that of the collapsed means and nobs is the right count.
  loglik <- -nobs * (opt$objective + cnst) - nobs * log(y_scale)

  # Extended BIC (Chen and Chen, 2008): bic plus 2*xi*k*log(p), where k is the
  # number of active jump locations and p = n-1 the number of candidates.
  # Corrects for searching over many candidate locations under sparsity --
  # plain bic's log(n)*df penalty does not charge for *which* locations were
  # selected, only for the resulting smoothing flexibility.
  k_active <- sum(opt$solution[2:n] > quasizero)
  bic_val  <- df * log(n) - 2 * loglik
  ebic_val <- bic_val + 2 * ebic_xi * k_active * log(n - 1)

  res <- list(
    opt                = opt,
    nobs               = n,
    df                 = df,
    maxsum             = maxsum_orig,
    loglik             = loglik,
    pars               = c(sigma  = opt$solution[1L] * y_scale / x_range^1.5,
                           lambda = lambda_orig),
    gamma              = opt$solution[2L:n] * x_range,
    ic                 = c(aic  =  2 * (df - loglik),
                           aicc =  2 * (df * n / (n - df - 1) - loglik),
                           bic  =  bic_val,
                           hq   =  2 * (log(log(n)) * df - loglik),
                           ebic =  ebic_val),
    smoothed_level     = (a1 + p11 * r1 + p12 * r2)[-(n + 1L)] * y_scale,
    smoothed_slope     = (a2 + p12 * r1 + p22 * r2)[-(n + 1L)] *
                           y_scale / x_range,
    var_smoothed_level = (p11 - p11^2 * n11 - 2 * p11 * p12 * n12 -
                            p12^2 * n22)[-(n + 1L)] * y_scale^2,
    x                  = x * x_range + x_min_s
  )
  new_ssj(res, method = "kfs", y = y * y_scale, call = match.call())
}


#' Automatic selection of the optimal time-warping smoothing spline with jumps
#' using the Kalman filter and smoother
#'
#' Runs \code{\link{ssj_kfs}} over a grid of \code{maxsum} values and returns
#' the solution with the lowest value of the selected information criterion.
#' Each grid point is warm-started from the previous optimum, and the search
#' stops early when the criterion increases for three consecutive grid values.
#'
#' @param y numeric vector with the response (missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param lambda positive scalar with the smoothing parameter.
#' @param grid numeric vector of \code{maxsum} values to search; the default is
#'   21 equally spaced values from 0 to the range of \code{x}.
#' @param ic string with the information criterion used for the choice: the
#'   default is \code{"ebic"} (recommended when jumps are expected to be rare
#'   relative to the \code{n-1} candidate locations, see \code{ebic_xi}), but
#'   \code{"bic"}, \code{"hq"}, \code{"aic"} and \code{"aicc"} are also
#'   available.
#' @param weights optional numeric vector of positive observation weights,
#'   passed to the fitter; see its help page.
#' @param edf logical: passed to \code{\link{ssj_kfs}} (default \code{TRUE}).
#' @param last_delta numeric scalar: passed to \code{\link{ssj_kfs}}
#'   (default 1).
#' @param ebic_xi numeric scalar in \code{[0,1]} passed to
#'   \code{\link{ssj_kfs}}: strength of the extended-BIC search penalty
#'   (default 0.15; only relevant when \code{ic = "ebic"}). See the
#'   "Choosing ebic_xi" section of \code{\link{ssj_kfs}} for how it was
#'   calibrated, and why it differs from the value used by
#'   \code{\link{auto_ssj_mle}}.
#'
#' @returns The output of the \code{\link{ssj_kfs}} function corresponding to
#'   the best choice according to the selected information criterion.
#'
#' @seealso \code{\link{ssj_kfs}}.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_kfs <- auto_ssj_kfs(y = faithful$waiting,
#'                        x = faithful$eruptions,
#'                        lambda = 1)
#' lines(x = of_kfs$x, y = of_kfs$smoothed_level, col = "blue", lwd = 2)
#'
#' @export
auto_ssj_kfs <- function(y, x, lambda, weights = NULL,
                         grid = seq(0, max(x) - min(x), length.out = 21L),
                         ic = c("ebic", "bic", "hq", "aic", "aicc"),
                         edf = TRUE, last_delta = 1, ebic_xi = 0.15) {
  ic          <- match.arg(ic)
  last_ic     <- Inf
  n_increases <- 0L
  parinit     <- NULL
  for (M in grid) {
    out <- ssj_kfs(y = y, x = x, lambda = lambda, weights = weights,
                   maxsum = M, edf = edf, parinit = parinit,
                   last_delta = last_delta, ebic_xi = ebic_xi)
    parinit    <- c(out$pars["sigma"], out$gamma)   # warm start
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
  best$call <- match.call()
  best
}
