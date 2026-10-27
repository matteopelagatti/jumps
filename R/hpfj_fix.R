#' HP filter with automatic jump detection and fixed smoothing constant
#'
#' Lower-level function for the Hodrick-Prescott filter with jumps and a fixed
#' smoothing parameter. Most users should call \code{\link{hpj}} instead, unless
#' they need more control and speed. Jumps happen contextually in the level and
#' in the slope: the standard deviation of the slope disturbance is
#' \eqn{\gamma} times the standard deviation of the level disturbance at time
#' \eqn{t}. Here the HP smoothing parameter \eqn{\lambda} is fixed by the user,
#' so that \eqn{\sigma^2_\varepsilon = \lambda\sigma^2_\zeta}.
#'
#' @param y numeric vector with the time series (missing values allowed).
#' @param lambda either a numeric scalar with the smoothing constant, or a
#'   string with the frequency of the time series, one of \code{"daily"},
#'   \code{"weekly"}, \code{"monthly"}, \code{"quarterly"}, \code{"annual"}; in
#'   the latter case the smoothing constant is computed according to Ravn and
#'   Uhlig (2002), that is \eqn{6.25 s^4}, where \eqn{s} is the number of
#'   observations per year.
#' @param maxsum maximum sum of additional level standard deviations.
#' @param edf logical: if \code{TRUE} (default) the effective degrees of freedom
#'   are computed, otherwise the number of degrees of freedom is computed in the
#'   LASSO-regression way.
#' @param parinit either \code{NULL} or a numeric vector of \code{2+n} starting
#'   values for the optimizer; the order of the parameters is sd(slope
#'   disturbance), square root of gamma, and \code{n} additional standard
#'   deviations for the level.
#'
#' @returns A list with the following slots:
#' \itemize{
#'  \item \code{opt}: the output of the optimization function (\code{nloptr}).
#'  \item \code{nobs}: number of observations.
#'  \item \code{df}: model's degrees of freedom.
#'  \item \code{maxsum}: the \code{maxsum} argument.
#'  \item \code{loglik}: value of the log-likelihood at the maximum.
#'  \item \code{pars}: named vector (sigma_slope, sigma_noise, gamma).
#'  \item \code{sigmas}: additional standard deviations of the level
#'    disturbance.
#'  \item \code{weights}: effective-degrees-of-freedom weights, or \code{NULL}
#'    when \code{edf = FALSE}.
#'  \item \code{ic}: vector of information criteria (aic, aicc, bic, hq).
#'  \item \code{smoothed_level}: vector with the smoothed level with jumps (HP
#'    filter with jumps).
#'  \item \code{var_smoothed_level}: variance of the smoothed level.
#' }
#'
#' @references
#' Maranzano, P. and Pelagatti, M. (2025). A Hodrick-Prescott filter with
#' automatically selected breaks. \emph{Economic Modelling}, \strong{150},
#' 107132. \doi{10.1016/j.econmod.2025.107132}
#'
#' Ravn, M. O. and Uhlig, H. (2002). On adjusting the Hodrick-Prescott filter
#' for the frequency of observations. \emph{The Review of Economics and
#' Statistics}, \strong{84}(2), 371--376. \doi{10.1162/003465302317411604}
#'
#' @seealso \code{\link{hpj}}, \code{\link{auto_hpfj_fix}}, \code{\link{hpfj}}.
#'
#' @examples
#' set.seed(202311)
#' n <- 100
#' mu <- 100*cos(3*pi/n*(1:n)) - ((1:n) > 50)*n - c(rep(0, 50), 1:50)*10
#' y <- mu + rnorm(n, sd = 20)
#' plot(y, type = "l")
#' lines(mu, col = "blue")
#' hp <- hpfj_fix(y, 60, 60)
#' lines(hp$smoothed_level, col = "red")
#'
#' @export
hpfj_fix <- function(y, lambda, maxsum = sd(y, na.rm = TRUE), edf = TRUE,
                     parinit = NULL) {
  if (is.numeric(lambda)) lambda <- abs(lambda)
  if (is.character(lambda)) {
    lambda <- match.arg(lambda[1],
                        c("daily", "weekly", "monthly", "quarterly", "annual"))
    lambda <- switch(lambda,
                     daily     = 110930628906,
                     weekly    = 45697600,
                     monthly   = 129600,
                     quarterly = 1600,
                     annual    = 6.25,
                     NULL)
    if (is.null(lambda)) stop("no valid value for lambda")
  }
  y <- as.numeric(y)
  nobs <- sum(!is.na(y))
  n   <- length(y)
  vy  <- var(y, na.rm = TRUE)
  sdy <- sqrt(vy)
  v_eta  <- numeric(n)
  v_zeta <- numeric(n)
  c_eta_zeta <- numeric(n)
  v_eps  <- numeric(n)
  a1  <- rep(y[!is.na(y)][1], n+1)
  a2  <- rep(0.0, n+1)
  p11 <- rep(vy, n+1)
  p12 <- numeric(n+1)
  p22 <- rep(vy/10, n+1)
  k1  <- numeric(n)
  k2  <- numeric(n)
  i   <- numeric(n)
  f   <- numeric(n)
  r1  <- numeric(n+1)
  r2  <- numeric(n+1)
  n11 <- numeric(n+1)
  n12 <- numeric(n+1)
  n22 <- numeric(n+1)
  e   <- numeric(n)
  d   <- numeric(n)
  w   <- numeric(n)
  cnst <- log(2*pi)/2
  vt_eta_ndx <- 3:(n+2)
  
  ##### Object function to optimize (log-lik with gradients w.r.t. sd_eta sd_zeta sd_)
  obj <- function(pars, wgt = FALSE) {
    # Time-varying variance of the level
    v_eta[]  <- pars[vt_eta_ndx]*pars[vt_eta_ndx]
    # Time-varying variance of the slope with chain/link to the level's variance
    v_zeta[] <- pars[1]*pars[1] + pars[2]*pars[2]*v_eta
    # Time-fixed variance of the measurement errors
    v_eps[]  <- pars[1]*pars[1]*lambda
    if (wgt) { # compute weights for edf
      mloglik <- -llt(y, v_eps, v_eta, v_zeta, c_eta_zeta,
                      a1, a2, p11, p12, p22,
                      k1, k2, i, f, r1, r2,
                      n11, n12, n22, e, d, w)
    } else { # no weights computed
      mloglik <- -llt(y, v_eps, v_eta, v_zeta, c_eta_zeta,
                      a1, a2, p11, p12, p22,
                      k1, k2, i, f, r1, r2,
                      n11, n12, n22, e, d)
    }
    rn1 <- r1[-1]*r1[-1] - n11[-1]
    rn2 <- r2[-1]*r2[-1] - n22[-1]
    list(
      # Mean/average log-likelihood for computational stability
      objective = mloglik/nobs,
      # Average gradient
      gradient  = c(
        # w.r.t sigma_zeta
        -sum(rn2)*pars[1]-sum(e*e-d)*pars[1]*lambda,
        # w.r.t gamma
        -sum(rn2*pars[2]*v_eta),
        # w.r.t. sigma_eta_t
        -(rn1 + rn2*pars[2]*pars[2])*pars[vt_eta_ndx]
      )/nobs
    )
  }
  
  ##### Constraints to the object function
  g <- function(pars, wgt) {
    list(
      constraints = sum(pars[vt_eta_ndx]) - maxsum,
      jacobian = c(0, 0, rep(1, n))
    )
  }
  
  ##### Optimization step
  ## Starting values
  if (is.null(parinit)) {
    inits <- c(sd_zeta = sdy/10, sqrt_gamma = 1/10, rep(1, n-1), 0)
  } else {
    inits <- parinit
  }
  ## Check on starting values
  # lb <- c(0, 0, rep(0, n))
  quasizero <- sdy*1.0e-9
  lb <- c(quasizero, quasizero, rep(0, n))
  inits[inits < lb] <- lb[inits < lb]
  
  # return(c(obj(inits),
  #          list(nd = numDeriv::grad(function(x) obj(x)[[1]], inits))))
  ## Optimization with CCSA ("conservative convex separable approximation")
  # see https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/
  opt <- nloptr::nloptr(x0 = inits,
                        eval_f = obj,
                        lb = lb,
                        eval_g_ineq = g,
                        opts = list(algorithm = "NLOPT_LD_CCSAQ",
                                    xtol_rel = 1.0e-5,
                                    check_derivatives = FALSE,
                                    maxeval = 2000),
                        wgt = FALSE
  )
  if (edf == TRUE) {
    obj(opt$solution, TRUE)
  }
  if (edf == TRUE) {
    df <- sum(w)
  } else {
    df <- 2 + sum(opt$solution[vt_eta_ndx] > 0)
  }
  loglik <- -nobs*(opt$objective + cnst)
  
  ##### Output list
  list(opt = opt,
       nobs = n,
       df = df,
       maxsum = maxsum,
       loglik = loglik,
       # NB: with a fixed lambda the parameter vector is
       # (sd_zeta, sqrt_gamma, sigma_eta_1, ..., sigma_eta_n), so gamma is
       # the square of solution[2] (solution[3] is the first sigma_eta).
       pars = c(sigma_slope = opt$solution[1],
                sigma_noise = opt$solution[1]*sqrt(lambda),
                gamma = opt$solution[2]*opt$solution[2]),
       sigmas = c(NA, opt$solution[-c(1, 2, n+2)]),
       weights = if (edf) w else NULL,
       ic = c(aic  = 2*(df - loglik),
              aicc = 2*(df*n/(n-df-1) - loglik),
              # bic  = log(df)*n - 2*loglik,
              bic  = df*log(n) - 2*loglik,
              hq   = 2*(log(log(n))*df - loglik)),
       smoothed_level = (a1 + p11*r1 + p12*r2)[-(n+1)],
       var_smoothed_level = (p11 - p11*p11*n11 - 2*p11*p12*n12 - p12*p12*n22)[-(n+1)]
  )
}



#' Automatic selection of the optimal HP filter with jumps and fixed smoothing
#' constant
#'
#' The regularization constant for the HP filter with jumps is the maximal sum
#' of standard deviations of the level disturbance. This value has to be passed
#' to the \code{\link{hpfj_fix}} function. \code{auto_hpfj_fix} runs
#' \code{hpfj_fix} on a grid of regularization constants and returns the output
#' of \code{hpfj_fix} selected by the chosen information criterion.
#'
#' @param y numeric vector containing the time series (missing values allowed).
#' @param lambda smoothing constant; see \code{\link{hpfj_fix}} for the
#'   accepted values.
#' @param grid numeric vector containing the grid for the argument
#'   \code{maxsum} of the \code{\link{hpfj_fix}} function.
#' @param ic string with the information criterion used for the choice: the
#'   default is \code{"bic"} (simulations show this is the best choice), but
#'   \code{"hq"}, \code{"aic"} and \code{"aicc"} are also available.
#' @param edf logical: \code{TRUE} (default) if the number of degrees of freedom
#'   should be computed as "effective degrees of freedom" (Efron, 1986), as
#'   opposed to a more traditional way (although not supported by theory) when
#'   \code{FALSE}.
#'
#' @returns The output of the \code{\link{hpfj_fix}} function corresponding to
#'   the best choice according to the selected information criterion.
#'
#' @references
#' Maranzano, P. and Pelagatti, M. (2025). A Hodrick-Prescott filter with
#' automatically selected breaks. \emph{Economic Modelling}, \strong{150},
#' 107132. \doi{10.1016/j.econmod.2025.107132}
#'
#' @seealso \code{\link{hpfj_fix}}, \code{\link{hpj}}.
#'
#' @examples
#' mod <- auto_hpfj_fix(Nile, "annual")
#' plot(as.numeric(Nile), type = "l")
#' lines(mod$smoothed_level, col = "red")
#'
#' @export
auto_hpfj_fix <- function(y, lambda,
                          grid = seq(0, sd(y, na.rm = TRUE)*10,
                                     sd(y, na.rm = TRUE)/10),
                          ic = c("bic", "hq", "aic", "aicc"),
                          edf = TRUE) {
  ic <- match.arg(ic)
  last_ic <- Inf
  for (M in grid) {
    out <- hpfj_fix(y, lambda, maxsum = M, edf = edf)
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
  best
}