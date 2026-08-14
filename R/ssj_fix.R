#' Smoothing splines with automatic jump detection and fixed smoothing
#' parameter
#'
#' Lower-level function for smoothing splines with jumps and a fixed smoothing
#' parameter. The additional variance of the level disturbance at each sampling
#' site is estimated by constrained maximum likelihood, using the Kalman filter
#' and smoother on the state-space form of the cubic smoothing spline.
#'
#' @param y numeric vector with the response (missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param lambda numeric scalar with the smoothing parameter.
#' @param maxsum maximum sum of the additional standard deviations. Note that
#'   this parameterises the jumps through the disturbance variances, so, unlike
#'   the budget of \code{\link{ssj_mle}}, \code{\link{ssj_kfs}} and
#'   \code{\link{ssj_atw}}, it is \emph{not} expressed in the units of \code{x}.
#' @param edf logical: if \code{TRUE} (default) the effective degrees of freedom
#'   are computed, otherwise the number of degrees of freedom is computed in the
#'   LASSO-regression way.
#' @param parinit either \code{NULL} or a numeric vector of \code{1+n} starting
#'   values for the optimizer; the order of the parameters is the constant
#'   standard deviation followed by \code{n} additional standard deviations.
#' @param last_delta numeric scalar with the difference \eqn{x_{n+1}-x_{n}}; it
#'   is not relevant for most applications and defaults to 1.
#'
#' @returns An object of class \code{ssj} (a list, so that all the slots below
#'   remain directly accessible) with a \code{\link[=print.ssj]{print}} and a
#'   \code{\link[=plot.ssj]{plot}} method. Its slots are:
#' \itemize{
#'  \item \code{opt}: the output of the optimization function (\code{nloptr}).
#'  \item \code{nobs}: number of observations.
#'  \item \code{df}: model's degrees of freedom.
#'  \item \code{maxsum}: the \code{maxsum} argument.
#'  \item \code{loglik}: value of the log-likelihood at the maximum.
#'  \item \code{pars}: named vector (sigma, lambda).
#'  \item \code{sigmas}: additional standard deviations of the level
#'    disturbance.
#'  \item \code{weights}: effective-degrees-of-freedom weights, or \code{NULL}
#'    when \code{edf = FALSE}.
#'  \item \code{ic}: vector of information criteria (aic, aicc, bic, hq).
#'  \item \code{smoothed_level}: vector with the smoothing spline with jumps.
#'  \item \code{var_smoothed_level}: vector with the variances of the smoothing
#'    spline.
#'  \item \code{x}: the ordered \code{x} variable.
#' }
#'
#' @seealso \code{\link{auto_ssj_fix}}, \code{\link{ssj_kfs}},
#'   \code{\link{ssj_mle}}, \code{\link{ssj_atw}}.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- ssj_fix(y = faithful$waiting,
#'                   x = faithful$eruptions,
#'                   lambda = 1,
#'                   maxsum = 150)
#' lines(x = of_jum$x, y = of_jum$smoothed_level, col = "red", lwd = 2)
#'
#' @export
ssj_fix <- function(y, x, lambda, maxsum = sd(y, na.rm = TRUE)/mean(diff(x)),
                    edf = TRUE, parinit = NULL, last_delta = 1) {
  n <- length(y)
  n1 <- n+1
  nobs <- sum(!is.na(y))
  vy <- var(y, na.rm = TRUE)
  sdy <- sqrt(vy)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  # NB: computed AFTER sorting.  var(diff(y)) is not invariant to the order in
  # which the observations are supplied, and it scales the diffuse prior on the
  # slope (p22 below), so taking it on the unsorted y would make the fit depend
  # on the row order of the input.
  vdy <- var(diff(y), na.rm = TRUE)
  delta <- c(diff(x), last_delta)
  
  var_eps <- numeric(n)
  var_eta <- numeric(n)
  var_zeta <- numeric(n)
  cov_eta_zeta <- numeric(n)
  a1  <- rep(y[!is.na(y)][1], n1)
  a2  <- numeric(n1)
  p11 <- rep(vy*1.0e5, n1)
  p12 <- numeric(n1)
  p22 <- rep(vdy*1.0e5, n1)
  k1  <- numeric(n)
  k2  <- numeric(n)
  i   <- numeric(n)
  f   <- numeric(n)
  r1  <- numeric(n1)
  r2  <- numeric(n1)
  n11 <- numeric(n1)
  n12 <- numeric(n1)
  n22 <- numeric(n1)
  e   <- numeric(n1)
  d   <- numeric(n1)
  w   <- numeric(n)
  cnst <- log(2*pi)/2
  
  D <- array(0, c(2, 2, n))
  D[1, 1, ] <- delta^3/3
  D[1, 2, ] <- D[2, 1, ] <- delta^2/2
  D[2, 2, ] <- delta
  
  # pars[1] sigma
  # pars[2:(n+1)] sigma_t
  vt_eta_ndx <- 2:(n+1)
  obj <- function(pars, wgt = FALSE) {
    sig2  <- pars[1]*pars[1] # common variance
    sig2t <- pars[vt_eta_ndx]*pars[vt_eta_ndx] # additional variances
    sum_sig2t <- sig2 + sig2t # common variance + additional variances
    var_eps[]  <- lambda*sig2 # variance of the noise
    var_eta[]  <- D[1, 1, ]*sum_sig2t # variance of eta
    var_zeta[] <- D[2, 2, ]*sum_sig2t # variance of zeta
    cov_eta_zeta[] <- D[1, 2, ]*sum_sig2t # covariance of eta and zeta
    
    if (wgt) { # compute weights for edf
      mloglik <- -llt_delta(y, delta, var_eps, var_eta, var_zeta, cov_eta_zeta,
                            a1, a2, p11, p12, p22,
                            k1, k2, i, f, r1, r2,
                            n11, n12, n22, e, d, w)
    } else { # no weights computed
      mloglik <- -llt_delta(y, delta, var_eps, var_eta, var_zeta, cov_eta_zeta,
                            a1, a2, p11, p12, p22,
                            k1, k2, i, f, r1, r2,
                            n11, n12, n22, e, d)
    }
    # Koopman-Shephard (2003) analytical scores via the smoother output
    rn1  <- (r1[-1]*r1[-1] - n11[-1])
    rn2  <- (r2[-1]*r2[-1] - n22[-1])
    rn12 <- (r1[-1]*r2[-1] - n12[-1])
    rn_sum <- rn1*D[1, 1, ] + 2*rn12*D[1, 2, ] + rn2*D[2, 2, ]
    list(
      # Minus mean log-likelihood for computational stability
      objective = mloglik/nobs,
      # Average gradient
      gradient  = c(
        # w.r.t sigma_zeta
        -pars[1]*sum(rn_sum) - pars[1]*lambda*sum(e*e-d),
        # w.r.t. sigma_eta_t
        -pars[vt_eta_ndx]*rn_sum
      )/nobs
    )
    
  }
  
  ##### Constraints to the object function
  g <- function(pars, wgt) {
    list(
      constraints = sum(pars[vt_eta_ndx]) - maxsum, # constraint <= 0
      jacobian = c(0, rep(1, n)) # derivatives of the constraint w.r.t. sigma_t
    )
  }
  
  ##### Optimization step
  ## Starting values
  if (is.null(parinit)) {
    inits <- c(sd_zeta = sdy/10, rep(1, n-1), 0) #
  } else {
    inits <- parinit
  }
  ## Check on starting values
  # lb <- c(0.0000000001, rep(0, n))
  quasizero <- sdy*1.0e-9
  lb <- c(quasizero, rep(0, n)) # lower bound
  inits[inits < lb] <- lb[inits < lb]      # make init coherent with lower bound

  ## Optimization with CCSA ("conservative convex separable approximation")
  # see https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  opt <- nloptr::nloptr(x0 = inits,
                        eval_f = obj,
                        lb = lb,
                        eval_g_ineq = g, # constraint
                        opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                    xtol_rel = 1.0e-5,
                                    check_derivatives = FALSE,
                                    maxeval = 2000),
                        wgt = FALSE # passed on to obj() and g() through ...
  )
  if (edf == TRUE) {
    obj(opt$solution, TRUE) # update pointed variables using estimated values
    df <- sum(w)
  } else {
    df <- 1 + sum(opt$solution[vt_eta_ndx] > quasizero) # LASSO-like count
  }
  loglik <- -nobs*(opt$objective + cnst)
  
  ##### Output list
  out <- list(opt = opt,
       nobs = n,
       df = df,
       maxsum = maxsum,
       loglik = loglik,
       pars = c(sigma = opt$solution[1], lambda = lambda),
       # the first element is empty because the shocks enter with a lag
       sigmas = c(NA, opt$solution[-c(1, n+1)]),
       weights = if (edf) w else NULL,
       ic = c(aic  = 2*(df - loglik),
              aicc = 2*(df*n/(n-df-1) - loglik),
              bic  = df*log(n) - 2*loglik,
              hq   = 2*(log(log(n))*df - loglik)),
       smoothed_level = (a1 + p11*r1 + p12*r2)[-(n+1)],
       var_smoothed_level = (p11 - p11*p11*n11 - 2*p11*p12*n12 - p12*p12*n22)[-(n+1)],
       x = x
  )
  new_ssj(out, method = "fix", y = y, call = match.call())
}

#' Automatic selection of the optimal smoothing spline with jumps and fixed
#' smoothing parameter
#'
#' The regularization parameter is the maximal sum of standard deviations added
#' to the static standard deviation. This value has to be passed to the
#' \code{\link{ssj_fix}} function. \code{auto_ssj_fix} runs \code{ssj_fix} on a
#' grid of regularization constants and returns the output of \code{ssj_fix}
#' selected by the chosen information criterion.
#'
#' @param y numeric vector with the response (missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param lambda numeric scalar with the smoothing parameter.
#' @param grid numeric vector containing the grid for the argument
#'   \code{maxsum} of the \code{\link{ssj_fix}} function.
#' @param ic string with the information criterion used for the choice: the
#'   default is \code{"bic"} (simulations show this is the best choice), but
#'   \code{"hq"}, \code{"aic"} and \code{"aicc"} are also available.
#' @param edf logical: \code{TRUE} (default) if the number of degrees of freedom
#'   should be computed as "effective degrees of freedom" (Efron, 1986), as
#'   opposed to a more traditional way (although not supported by theory) when
#'   \code{FALSE}.
#' @param last_delta numeric scalar with the difference \eqn{x_{n+1}-x_{n}}; not
#'   relevant for most applications and defaults to 1.
#'
#' @returns The output of the \code{\link{ssj_fix}} function corresponding to
#'   the best choice according to the selected information criterion.
#'
#' @seealso \code{\link{ssj_fix}}.
#'
#' @examples
#' # a short, explicit grid keeps the example fast; the default grid has 101
#' # points and is correspondingly slower
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- auto_ssj_fix(y = faithful$waiting,
#'                        x = faithful$eruptions,
#'                        lambda = 1,
#'                        grid = seq(0, 400, length.out = 5))
#' lines(x = of_jum$x, y = of_jum$smoothed_level, col = "red", lwd = 2)
#'
#' @export
auto_ssj_fix <- function(y, x, lambda,
                         grid = seq(0, sd(y, na.rm = TRUE)/mean(diff(x))*10,
                                    sd(y, na.rm = TRUE)/mean(diff(x))/10),
                         ic = c("bic", "hq", "aic", "aicc"),
                         edf = TRUE, last_delta = 1) {
  ic <- match.arg(ic)
  last_ic <- Inf
  for (M in grid) {
    out <- ssj_fix(y = y, x = x, lambda = lambda,
                   maxsum = M, edf = edf, last_delta = last_delta)
    current_ic <- switch (ic,
                          bic = out$ic["bic"],
                          hq = out$ic["hq"],
                          aic = out$ic["aic"],
                          aicc = out$ic["aicc"]
    )
    if (current_ic < last_ic) {
      best <- out
      last_ic <- current_ic
    }
  }
  best$call <- match.call()
  best
}
