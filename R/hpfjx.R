#' HP filter with jumps and regressors (still experimental)
#'
#' This function needs more testing since it does not seem to work as expected.
#' For this reasin the wrapper \code{hpj} at the moment does not allow regressors.
#' This is the same as \code{hpfj} but with the possibility of including regressors.
#' The regressors should be zero-mean so that the HP filter can be interpreted as
#' a mean value of the time series. 
#' Jumps happen contextually in the level and in the slope: the standard deviation
#' of the slope disturbance is \eqn{\gamma} times the standard deviation of the
#' level disturbance at time \eqn{t}.
#' The HP smoothing parameter \eqn{\lambda} is estimated via MLE (assuming normally
#' distributed disturbances) as in Wahba (1978):
#' \eqn{\lambda = \sigma^2_\varepsilon / \sigma^2_\zeta}.
#' 
#' @references Whaba (1978)
#' "Improper priors, spline smoothing and the problem of guarding against model errors in regression",
#' *Journal of the Royal Statistical Society. Series B*, Vol. 40(3), pp. 364-372.
#' DOI:10.1111/j.2517-6161.1978.tb01050.x
#' 
#' @param y vector with the time series
#' @param X matrix with regressors in the columns
#' @param maxsum maximum sum of additional level standard deviations;
#' @param edf boolean if TRUE computes effective degrees of freedom otherwise computes
#' the number of degrees of freedom in the LASSO-regression way.
#' @param parinit either NULL or vector of 3+n parameters with starting values for the
#' optimizer; the order of the parameters is sd(slope disturbnce), sd(observatio noise),
#' square root of gamma, n additional std deviations for the slope
#' @return list with the following slots:
#' \itemize{
#'  \item opt: the output of the optimization function (nloptr)
#'  \item nobs: number of observations
#'  \item df: number of estimated parameters (model's degrees of freedom)
#'  \item loglik: value of the log-likelihood at maximum
#'  \item ic: vector of information criteria (aic, aicc, bic, hq)
#'  \item smoothed_level: vector with smoothed level with jumps (hp filter with jumps)
#'  \item var_smoothed_level: variance of the smoothed level
#' }
#' @examples 
#' y <- log(AirPassengers)
#' n <- length(y)
#' mod <- hpfjx(y, trigseas(n, 12))
#' hpj <- ts(mod$smoothed_level, start(y), frequency = 12)
#' plot(y)
#' lines(hpj, col = "red")
#' 
#' @export
hpfjx <- function(y, X, maxsum = sd(y), edf = TRUE, parinit = NULL) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  k <- ncol(X) # number of regressors
  n   <- length(y) # time series length
  vy  <- var(y, na.rm = TRUE)
  sdy <- sqrt(vy)
  v_eta  <- numeric(n)
  v_zeta <- numeric(n)
  c_eta_zeta <- numeric(n)
  v_eps  <- numeric(n)
  a1  <- rep(y[1], n+1)
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
  vt_eta_ndx <- 4:(n+3)
  vt_reg_ndx <- (n+4):(n+3+k)
  
  ##### Object function to optimize (minus mean log-lik with gradients)
  obj <- function(pars, wgt = FALSE) {
    # Time-varying variance of the level
    v_eta[]  <- pars[vt_eta_ndx]*pars[vt_eta_ndx]
    # Time-varying variance of the slope with chain/link to the level's variance
    v_zeta[] <- pars[1]*pars[1] + pars[3]*pars[3]*v_eta
    # Time-fixed variance of the measurement errors
    v_eps[]  <- pars[2]*pars[2]
    yx <- as.numeric(y - X %*% pars[vt_reg_ndx])
    if (wgt) { # compute weights for edf
      mloglik <- -llt(yx, v_eps, v_eta, v_zeta, c_eta_zeta,
                      a1, a2, p11, p12, p22,
                      k1, k2, i, f, r1, r2,
                      n11, n12, n22, e, d, w)
    } else { # no weights computed
      mloglik <- -llt(yx, v_eps, v_eta, v_zeta, c_eta_zeta,
                      a1, a2, p11, p12, p22,
                      k1, k2, i, f, r1, r2,
                      n11, n12, n22, e, d)
    }
    rn1 <- r1[-1]*r1[-1] - n11[-1]
    rn2 <- r2[-1]*r2[-1] - n22[-1]
    da1 <- matrix(0, n, k)
    da2 <- matrix(0, n, k)
    da(k1, k2, X, da1, da2)
    # the above Rcpp function does the following
    # for (t in 1:(n-1)) {
    #   da1[t+1, ] <- (1-k1[t])*da1[t, ] + da2[t, ] - k1[t]*X[t, ]
    #   da2[t+1, ] <- da2[t, ] - k2[t]*X[t, ] - k2[t]*da1[t, ]
    # }
    list(
      # Mean negative log-likelihood for computational stability
      objective = mloglik/n,
      # Average gradient
      gradient  = c(
        # w.r.t sigma_zeta
        -sum(rn2)*pars[1],
        # w.r.t sigma_epsilon
        -sum(e*e-d)*pars[2],
        # w.r.t gamma
        -sum(rn2*pars[3]*v_eta),
        # w.r.t. sigma_eta_t
        -(rn1 + rn2*pars[3]^2)*pars[vt_eta_ndx],
        # w.r.t. beta
        -crossprod(i/f, X+da1)
      )/n
    )
  }
  
  ##### Constraints to the object function
  g <- function(pars, wgt) {
    list(
      constraints = sum(pars[vt_eta_ndx]) - maxsum,
      jacobian = c(0, 0, 0, rep(1, n), rep(0, k))
    )
  }
  
  ##### Optimization step
  ## Starting values
  if (is.null(parinit)) {
    inits <- c(sd_zeta = sdy/10, sd_eps = sdy, sqrt_gamma = 1/10,
               rep(1, n-1), 0,
               rep(0, k))
  } else {
    inits <- parinit
  }
  ## Lower bounds
  # lb <- c(0, 0, 0, rep(0, n), rep(-Inf, k))
  quasizero <- sdy*1.0e-9
  lb <- c(quasizero, quasizero, quasizero, rep(0, n), rep(-Inf, k))
  nrndx <- 1:(n+3)
  inits[nrndx][inits[nrndx] < lb[nrndx]] <- lb[nrndx][inits[nrndx] < lb[nrndx]]
  
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
                                    maxeval = 500),
                        wgt = FALSE
  )
  if (edf == TRUE) {
    obj(opt$solution, TRUE)
    df <- sum(w) + k
  } else {
    df <- 3 + sum(opt$solution[vt_eta_ndx] > 0) + k
  }
  loglik <- -n*(opt$objective + cnst)
  
  ##### Output list
  list(opt = opt,
       nobs = n,
       df = df,
       maxsum = maxsum,
       loglik = loglik,
       coefficients = opt$solution[vt_reg_ndx],
       pars = c(sigma_slope = opt$solution[1],
                sigma_noise = opt$solution[2],
                gamma = opt$solution[3]*opt$solution[3]),
       sigmas = c(NA, opt$solution[vt_eta_ndx][-n]),
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

#' Automatic selection of the optimal HP filter with jumps and regressors
#' 
#' This function needs more testing since it does not seem to work as expected.
#' For this reasin the wrapper \code{hpj} at the moment does not allow regressors.
#' The regularization constant for the HP filter with jumps is the
#' maximal sum of standard deviations for the level disturbance. This value
#' has to be passed to the \code{hpfjx} function. The \code{auto_hpfjx} runs
#' \code{hpfjx} on a grid of regularization constants and returns the output
#' of \code{hpfjx} selected by the chosen information criterion.
#' 
#' @param y numeric vector cotaining the time series;
#' @param X numeric matrix with regressors in the columns;
#' @param grid numeric vector containing the grid for the argument \code{maxsum}
#' of the \code{hpfj} function;
#' @param ic string with information criterion for the choice: the default is
#' "bic" (simulations show this is the best choice), but also "hq", "aic" and "aicc"
#' are available;
#' @param edf logical scalar: TRUE (default) if the number of degrees of freedom
#' should be computed as "effective degrees of freedom" (Efron, 1986) as opposed
#' to a more traditional way (although not supported by theory) when FALSE.
#' 
#' @returns The ouput of the \code{hpjf} function corresponding to the best
#' choice according to the selected information criterion.
#' 
#' @examples 
#' y <- log(AirPassengers)
#' n <- length(y)
#' mod <- auto_hpfjx(y, trigseas(n, 12))
#' hpj <- ts(mod$smoothed_level, start(y), frequency = 12)
#' plot(y)
#' lines(hpj, col = "red")
#' 
#' @export
auto_hpfjx <- function(y, X, grid = seq(0, sd(y)*10, sd(y)/10),
                      ic = c("bic", "hq", "aic", "aicc"), edf = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (length(y) != dim(X)[1]) stop("y and X have a different number of observations")
  ic <- match.arg(ic)
  best_ic <- Inf
  for (M in grid) {
    out <- hpfjx(y = y, X = X, maxsum = M, edf = edf)
    current_ic <- switch (ic,
                          bic  = out$ic["bic"],
                          hq   = out$ic["hq"],
                          aic  = out$ic["aic"],
                          aicc = out$ic["aicc"]
    )
    if (current_ic < best_ic) {
      best <- out
      best_ic <- current_ic
    }
  }
  best
}
