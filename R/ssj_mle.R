#' Smoothing splines with automatic jump detection via time warping and
#' maximum-likelihood estimation of the smoothing parameter
#'
#' Identical to \code{\link{ssj_kfs}}, except that \code{lambda} is no longer a
#' fixed input: it is a free parameter optimized jointly with \eqn{\sigma} and
#' \eqn{\gamma}. This avoids having to grid-search \code{lambda} (which would
#' suffer from the same degrees-of-freedom problem as \code{maxsum}) on top of
#' the \code{maxsum} grid search.
#'
#' Since \code{lambda} enters only the observation variance \eqn{H_t =
#' \lambda\sigma^2} and never the transition matrix, its score does not need the
#' extra correction that \eqn{\gamma} requires: the plain Koopman-Shephard
#' formula applies directly, \eqn{\partial\ell/\partial\lambda =
#' (\sigma^2/2)\sum_t(e_t^2-d_t)}. Internally \code{lambda} is optimized via
#' \eqn{\mu=\log\lambda}, which guarantees \eqn{\lambda > 0} without a box
#' constraint and avoids the scale instability observed when optimizing the raw
#' \code{lambda} directly (its sensible range can span orders of magnitude, as
#' for any variance ratio); the score is then a one-line chain-rule correction,
#' \eqn{\partial\ell/\partial\mu = \lambda\,\partial\ell/\partial\lambda}. Both
#' \code{parinit} and the returned \code{pars} are unaffected and stay in terms
#' of the raw \code{lambda}.
#'
#' When \code{parinit} is not supplied, the success of the optimizer depends on
#' the path it traverses, not just on how close the starting \code{lambda} is to
#' the eventual optimum. Rather than relying on a single formula,
#' \code{ssj_mle} therefore tries a fixed grid of starting values for
#' \code{lambda} spanning several orders of magnitude and keeps whichever
#' converges to the highest likelihood.
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
#' @param maxsum maximum sum of the \eqn{\gamma_t}, that is the warping budget
#'   \eqn{M} of the constraint \eqn{\sum_t \gamma_t \le M}. It is expressed in
#'   the units of \code{x} and defaults to the range of \code{x}, the same
#'   upper end as the grid searched by \code{\link{auto_ssj_mle}}.
#' @param weights optional numeric vector of positive observation weights, one
#'   per element of \code{y}. The observation variance becomes
#'   \eqn{H_t = \lambda\sigma^2/w_t}, so a point carrying weight \eqn{w_t} is
#'   treated as the mean of \eqn{w_t} independent observations. This is what
#'   makes it correct to collapse tied values of \code{x} to their mean and
#'   pass the number of replicates as the weight: the Gaussian likelihood
#'   factorises, so the weighted fit on the collapsed data is identical to the
#'   fit on the raw replicates. The default, \code{NULL}, means unit weights.
#'
#'   Note a caveat specific to this function, which estimates the noise level
#'   rather than taking it as given. Collapsing replicates to their mean
#'   discards the within-group scatter, and that scatter is exactly what
#'   identifies the observation variance: on a test case with 150 observations
#'   at 40 distinct sites the raw fit recovered a variance of 0.021 against a
#'   true 0.0225, while the weighted collapsed fit returned 0.086. Pass the
#'   discarded scatter through \code{ss_within} to restore it; the estimates
#'   then agree with the raw ones to four or five significant digits. At a
#'   \emph{fixed} smoothing parameter the discarded term does not involve
#'   \eqn{f} at all, which is why \code{\link{ssj_kfs}},
#'   \code{\link{ssj_fix}} and \code{\link{ssj_atw}} reproduce the raw fit
#'   without it.
#' @param ss_within optional pooled within-group sum of squares discarded when
#'   replicates were collapsed to their means, that is
#'   \eqn{\sum_i\sum_j (y_{ij}-\bar y_i)^2}. Supplying it, together with the
#'   replicate counts in \code{weights}, restores the term of the likelihood
#'   that the collapsing removed, so the estimates are those that the raw
#'   replicates would have produced. It has no effect on the fit at a fixed
#'   smoothing parameter; what it restores is the information identifying the
#'   observation variance, and hence \code{lambda}.
#' @param edf logical: if \code{TRUE} (default) the effective degrees of freedom
#'   are computed as the trace of the hat matrix; if \code{FALSE} the non-zero
#'   \eqn{\gamma_t} are counted.
#' @param parinit \code{NULL} or a numeric vector of length \code{n+1} with
#'   starting values for the optimizer, namely \eqn{(\sigma, \lambda, \gamma_1,
#'   \ldots, \gamma_{n-1})}; pass the previous solution for warm starting.
#' @param last_delta positive scalar: dummy spacing appended after the last
#'   observation (it does not affect the fit; default 1).
#' @param ebic_xi numeric scalar in \code{[0,1]}: strength of the extended-BIC
#'   penalty for searching over the \code{n-1} candidate jump locations
#'   (default 0.3; 0 reduces \code{ebic} to plain \code{bic}). The default is
#'   calibrated rather than taken from theory; see the note below.
#'
#' @section Choosing \code{ebic_xi}:
#' The default is the value that minimises worst-case regret in root mean
#' squared error over a factorial design of 48 cells -- four sample sizes,
#' three noise levels and four counts of discontinuities, 17,400 replications
#' in all. The experiment is shipped as \code{inst/scripts/xi_calibration.R}.
#' For this estimator the optimum is flat in \eqn{n}, at 0.27 to 0.33, so a
#' constant is appropriate; \code{\link{ssj_kfs}} and \code{\link{ssj_atw}}
#' were calibrated separately and carry different defaults, because they are
#' different procedures even where they fit the same model.
#'
#' The value sits below the Chen and Chen (2008) consistency floor of about
#' 0.5 for \eqn{p = n-1} candidates, and that is not an oversight. Their
#' derivation assumes the base penalty charges only for dimension, with the
#' whole cost of searching the model space carried by the \eqn{\xi} term. Here
#' part of it is already paid: \code{df} is the trace of the smoothing matrix
#' and grows with each activated warp, so the plain BIC term levies about
#' \eqn{0.16\log n} per declared jump before the extended term adds anything.
#' Expressed on the \eqn{\xi} scale that implicit charge is roughly 0.08.
#'
#' Detection accuracy cannot be used to calibrate \eqn{\xi} on this problem:
#' with \eqn{n-1} candidates and a handful of true discontinuities the
#' negative class dominates, and both balanced accuracy and the Matthews
#' coefficient give optima that wander across the design. Larger values
#' declare fewer jumps if that is what a particular application needs.
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
#'  \item \code{pars}: named vector (sigma, lambda); \code{lambda} is estimated.
#'  \item \code{gamma}: estimated time-warp vector of length \code{n-1}.
#'  \item \code{ic}: vector of information criteria (aic, aicc, bic, hq, ebic).
#'  \item \code{smoothed_level}: smoothed estimate of \eqn{f} at each \code{x}.
#'  \item \code{var_smoothed_level}: pointwise variance of the smoothed level.
#'  \item \code{dist_level}, \code{dist_slope}: smoothed level and slope
#'    disturbances (\code{NULL} when \code{edf = FALSE}).
#'  \item \code{std_dist_level}, \code{std_dist_slope}: the corresponding
#'    standardised (auxiliary residual) versions.
#'  \item \code{x}: the ordered \code{x} variable.
#' }
#'
#' @references
#' Pelagatti, M., Maranzano, P. and Gatto Ronchieri, G. C. (2026). Smoothing
#' splines with some discontinuities. \emph{arXiv preprint} arXiv:2610.33137.
#' \doi{10.48550/arXiv.2610.33137}
#'
#' Chen, J. and Chen, Z. (2008). Extended Bayesian information criteria for
#' model selection with large model spaces. \emph{Biometrika}, \strong{95}(3),
#' 759--771. \doi{10.1093/biomet/asn034}
#'
#' Koopman, S. J. and Shephard, N. (1992). Exact score for time series models in
#' state space form. \emph{Biometrika}, \strong{79}(4), 823--826.
#' \doi{10.1093/biomet/79.4.823}
#'
#' @seealso \code{\link{auto_ssj_mle}}, \code{\link{ssj_kfs}},
#'   \code{\link{ssj_atw}}, \code{\link{ssj_fix}}.
#'
#' @examples
#' # 'maxsum' is a warping budget in the units of x; the range of eruptions
#' # is about 3.5, so a budget of 1 is already generous.
#' # thinned for speed: the full 'faithful' data set works just as well
#' fa <- faithful[seq(1, nrow(faithful), by = 3), ]
#' plot(fa$eruptions, fa$waiting)
#' of_mle <- ssj_mle(y = fa$waiting, x = fa$eruptions, maxsum = 1)
#' lines(x = of_mle$x, y = of_mle$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
ssj_mle <- function(y, x, maxsum = max(x) - min(x), weights = NULL,
                    ss_within = NULL, edf = TRUE, parinit = NULL,
                    last_delta = 1, ebic_xi = 0.3) {
  wt <- ssj_check_weights(weights, length(y))
  y  <- as.numeric(y); x <- as.numeric(x)
  # Within-group scatter discarded when replicates were collapsed to their
  # mean.  It carries no information about f, but it is what identifies the
  # observation variance, so it must re-enter the likelihood for lambda to be
  # estimated as if the raw replicates had been supplied.
  if (is.null(ss_within)) {
    ssw <- 0; df_within <- 0
  } else {
    if (!is.numeric(ss_within) || length(ss_within) != 1L ||
        !is.finite(ss_within) || ss_within < 0)
      stop("'ss_within' must be a single non-negative number")
    df_within <- sum(wt) - length(y)
    if (df_within <= 0)
      stop("'ss_within' needs 'weights' summing to more than length(y): the ",
           "weights are the replicate counts behind the collapsed means")
    ssw <- ss_within
  }
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
  # add the scatter discarded by the automatic collapsing above
  ssw       <- ssw + cl_$ss_within
  df_within <- df_within + cl_$df_within
  ssw <- ssw / y_scale^2   # same units as the normalised y
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

    # Contribution of the discarded within-group scatter.  Writing
    # s2e = lambda*sigma^2 for the per-observation noise variance, the raw
    # log-likelihood exceeds the collapsed one by
    #   -0.5 * (df_within * log(s2e) + ss_within / s2e)
    # up to constants free of the parameters.
    if (df_within > 0) {
      s2e     <- lambda * sig2
      corr    <- -0.5 * (df_within * log(s2e) + ssw / s2e)
      mloglik <- mloglik - corr
      dcorr   <- -0.5 * (df_within / s2e - ssw / s2e^2)   # d corr / d s2e
    } else {
      dcorr <- 0
    }

    # Koopman-Shephard (2003) analytical scores via smoother output.
    rn1  <- r1[-1]^2       - n11[-1]
    rn2  <- r2[-1]^2       - n22[-1]
    rn12 <- r1[-1] * r2[-1] - n12[-1]

    # H_t = lambda*sigma^2 enters ONLY the observation variance, never the
    # transition matrix, so both sigma's and lambda's scores follow the plain
    # Koopman-Shephard formula with no extra correction term -- unlike gamma.
    process_s   <- sum(rn1 * tauf^3/3 + rn12 * tauf^2 + rn2 * tauf)
    obs_s       <- sum((e * e - d) / wf)          # raw, no lambda factor
    grad_sigma  <- sigma * (process_s + lambda * obs_s)   # dH_t/d(sigma)  = 2*lambda*sigma
    grad_lambda <- (sig2 / 2) * obs_s                      # dH_t/d(lambda) = sigma^2
    # chain rule for the within-group term: s2e = lambda*sigma^2, so
    # d s2e/d sigma = 2*lambda*sigma and d s2e/d lambda = sigma^2
    grad_sigma  <- grad_sigma  + dcorr * 2 * lambda * sigma
    grad_lambda <- grad_lambda + dcorr * sig2
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
  ssj_check_parinit(parinit, n + 1L,
                    "(sigma, lambda, gamma_1, ..., gamma_{n-1})",
                    n, cl_$collapsed)
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
  if (df_within > 0)
    loglik <- loglik - 0.5 * (df_within * log(2 * pi) + sum(log(wt)))
  # Back to the scale of the data as supplied.  y was divided by y_scale on
  # entry, and the density of the rescaled copy differs from the density of the
  # original by the Jacobian of that map: one factor of y_scale per observation
  # entering the likelihood.  Those are the nobs collapsed means plus, when the
  # within-group scatter is restored, its df_within degrees of freedom -- which
  # together are exactly the number of raw observations.
  loglik <- loglik - (nobs + df_within) * log(y_scale)

  # Extended BIC (Chen and Chen, 2008): bic plus 2*xi*k*log(p), where k is the
  # number of active jump locations and p = n-1 the number of candidates.
  k_active <- sum(opt$solution[3:(n + 1)] > quasizero)
  bic_val  <- df * log(n) - 2 * loglik
  ebic_val <- bic_val + 2 * ebic_xi * k_active * log(n - 1)

  # Back-transform all returned quantities to original units.
  res <- list(
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
    smoothed_slope     = (a2 + p12 * r1 + p22 * r2)[-(n + 1L)] *
                           y_scale / x_range,
    var_smoothed_level = (p11 - p11^2 * n11 - 2 * p11 * p12 * n12 -
                            p12^2 * n22)[-(n + 1L)] * y_scale^2,
    dist_level         = eta_hat  * y_scale,          # level disturbance in orig y units
    dist_slope         = zeta_hat * y_scale / x_range, # slope disturbance in orig y/x units
    std_dist_level     = std_eta_hat,                  # dimensionless, no rescaling
    std_dist_slope     = std_zeta_hat,
    x                  = x * x_range + x_min_s        # back to original x
  )
  new_ssj(res, method = "mle", y = y * y_scale, call = match.call())
}


#' Automatic selection of the optimal time-warping smoothing spline with jumps
#' using maximum-likelihood estimation
#'
#' Runs \code{\link{ssj_mle}} over a grid of \code{maxsum} values and returns
#' the solution with the lowest value of the selected information criterion.
#' Unlike \code{\link{auto_ssj_kfs}} and \code{\link{auto_ssj_atw}},
#' \code{lambda} is not part of the grid: it is estimated by maximum likelihood
#' inside every \code{ssj_mle} call, so that only \code{maxsum} has to be
#' searched.
#'
#' The pure smooth (\code{maxsum = 0}) is fitted first: its smoothed level
#' disturbances are large exactly where the smooth model is "surprised" by a
#' sudden level change, and are used to allocate the \code{maxsum} budget
#' proportionally, giving the optimizer a warm start at the right locations.
#' The search stops early when the criterion increases for three consecutive
#' grid values.
#'
#' @param y numeric vector with the response (missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param grid numeric vector of \code{maxsum} values to search; the default is
#'   21 equally spaced values from 0 to the range of \code{x}.
#' @param ic string with the information criterion used for the choice: the
#'   default is \code{"ebic"} (recommended when jumps are expected to be rare
#'   relative to the \code{n-1} candidate locations, see \code{ebic_xi}), but
#'   \code{"bic"}, \code{"hq"}, \code{"aic"} and \code{"aicc"} are also
#'   available.
#' @param weights optional numeric vector of positive observation weights,
#'   passed to the fitter; see its help page.
#' @param ss_within optional pooled within-group sum of squares, passed to
#'   \code{\link{ssj_mle}}; see its help page.
#' @param edf logical: passed to \code{\link{ssj_mle}} (default \code{TRUE}).
#' @param last_delta numeric scalar: passed to \code{\link{ssj_mle}}
#'   (default 1).
#' @param ebic_xi numeric scalar in \code{[0,1]} passed to
#'   \code{\link{ssj_mle}}: strength of the extended-BIC search penalty
#'   (default 0.3; only relevant when \code{ic = "ebic"}). See the
#'   "Choosing ebic_xi" section of \code{\link{ssj_mle}} for how it was
#'   calibrated.
#'
#' @returns The output of the \code{\link{ssj_mle}} function corresponding to
#'   the best choice according to the selected information criterion.
#'
#' @references
#' Pelagatti, M., Maranzano, P. and Gatto Ronchieri, G. C. (2026). Smoothing
#' splines with some discontinuities. \emph{arXiv preprint} arXiv:2610.33137.
#' \doi{10.48550/arXiv.2610.33137}
#'
#' @seealso \code{\link{ssj_mle}}.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_mle <- auto_ssj_mle(y = faithful$waiting,
#'                        x = faithful$eruptions)
#' lines(x = of_mle$x, y = of_mle$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
auto_ssj_mle <- function(y, x, weights = NULL, ss_within = NULL,
                         grid = seq(0, max(x) - min(x), length.out = 21L),
                         ic = c("ebic", "bic", "hq", "aic", "aicc"),
                         edf = TRUE, last_delta = 1, ebic_xi = 0.3) {
  ic <- match.arg(ic)

  # Stage 1: fit the pure smooth (maxsum = 0).  With no gamma degrees of
  # freedom the optimisation landscape is unimodal and sigma/lambda are
  # identified cleanly.  The smoothed level disturbances E[eta_t | Y] are
  # large at locations where the smooth model is "surprised" by a sudden
  # level change — i.e. at genuine discontinuities.
  fit0   <- ssj_mle(y = y, x = x, weights = weights, ss_within = ss_within,
                    maxsum = 0, edf = TRUE,
                    last_delta = last_delta, ebic_xi = ebic_xi)
  # ssj_mle() collapses tied values of x, so the fit is returned on the
  # DISTINCT sites: everything sized per observation below must follow that
  # length, not the length of the raw input.
  n      <- length(fit0$x)
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
                     weights = weights, ss_within = ss_within,
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
  best$call <- match.call()
  best
}
