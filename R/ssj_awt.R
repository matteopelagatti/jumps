#' Smoothing splines with with automatic jumps detection and
#' fixed smoothing parameter using the Alternating time-warping algorithm
#' 
#' Similar to ssj_fix, but using the Alternating time-warping algorithm. 
#' 
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param lambda a numeric scalar with the smoothing parameter.
#' @param maxsum maximum sum of additional standard deviations.
#' @param max_iter maximum number of iterations (default 100).
#' @param tol tolerance for convergence (default 1e-4).
#' @param learning_rate initial and fallback step size for the Barzilai-Borwein
#'   adaptive scheme (default 0.1).
#' @param gamma_init optional numeric vector of length n-1 with a warm-start
#'   for gamma; ignored when its length does not match n-1.
#' @param ebic_xi numeric scalar in [0,1]: strength of the Extended-BIC penalty
#'   for searching over the n-1 candidate jump locations (default 1, the
#'   Chen and Chen 2008 recommendation for sparse selection; 0 reduces ebic to
#'   plain bic).
#' @return list with the following slots:
#' \itemize{
#'  \item nobs: number of observations.
#'  \item df: model's degrees of freedom.
#'  \item maxsum: maximum sum of additional standard deviations.
#'  \item loglik: Gaussian log-likelihood evaluated at the unbiased residual variance sigma^2 = rss/(n - edf).
#'  \item gamma: additional times.
#'  \item ic: vector of information criteria (aic, aicc, bic, hq, ebic).
#'  \item smoothed_level: vector with smoothing splines with jumps.
#'  \item x: the ordered x variable
#' }
#' 
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- ssj_atw(y = faithful$waiting,
#'                   x = faithful$eruptions,
#'                   lambda = 1,
#'                   maxsum = 150)
#' lines(x = sort(faithful$eruptions),
#'       y = of_jum$smoothed_level, col = "red", lwd = 2)
#'       
#' @export
ssj_atw <- function(y, x, lambda, maxsum = sd(y, na.rm = TRUE)/mean(diff(x)),
                    max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                    gamma_init = numeric(0), ebic_xi = 1) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  ord <- order(x)
  y <- y[ord]
  x <- x[ord]
  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    y   <- as.numeric(tapply(y, grp, mean, na.rm = TRUE))
    x   <- x[!duplicated(x)]
  }
  if (length(gamma_init) != length(x) - 1L) gamma_init <- numeric(0)
  out <- solve_jump_spline_fast(x_in = x, y_in = y,
                                lambda = lambda, M = maxsum,
                                max_iter = max_iter, tol = tol,
                                learning_rate = learning_rate,
                                gamma_init = gamma_init,
                                ebic_xi = ebic_xi)
  list(nobs = length(y),
       df = out$edf,
       maxsum = maxsum,
       loglik = out$loglik,
       gamma = out$gamma,
       ic = c(aic  = out$aic,
              aicc = out$aicc,
              bic  = out$bic,
              hq   = out$hq,
              ebic = out$ebic),
       smoothed_level = out$f,
       x = x
  )
}

#' Automatic selection of the optimal smoothing spline function with jumps
#' and fixed smoothing parameter
#' 
#' The regularization parameter for the HP filter with jumps is the
#' maximal sum of additional times. This value
#' has to be passed to the \code{ssj_atw} function. The function \code{auto_ssj_atw}
#' runs \code{ssj_atw} on a grid of regularization constants and returns the
#' output of \code{ssj_atw} according to the chosen information criterion.
#' Each grid point is warm-started from the solution at the previous point.
#' The grid search stops early when the selected IC increases for three
#' consecutive grid values.
#' 
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param lambda smoothing constant;
#' @param grid numeric vector containing the grid for the argument \code{maxsum}
#' of the \code{ssj_atw} function;
#' @param ic string with information criterion for the choice: the default is
#' "bic", but "ebic" (recommended when jumps are expected to be rare relative
#' to the n-1 candidate locations -- see \code{ebic_xi}), "hq", "aic" and
#' "aicc" are also available;
#' @param max_iter maximum number of iterations in ssj_atw (default 100).
#' @param tol tolerance for convergence of ssj_atw (default 1e-4).
#' @param learning_rate learning rate of ssj_atw (default 0.1).
#' @param ebic_xi numeric scalar in [0,1] passed to \code{ssj_atw}: strength
#' of the Extended-BIC search penalty (default 1; only relevant when
#' \code{ic = "ebic"}).
#'
#' @returns The ouput of the \code{ssj_awt} function corresponding to the best
#' choice according to the selected information criterion.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- auto_ssj_atw(y = faithful$waiting,
#'                        x = faithful$eruptions,
#'                        lambda = 1)
#' lines(x = sort(faithful$eruptions),
#'       y = of_jum$smoothed_level, col = "red", lwd = 2)
#'
#' @export
auto_ssj_atw <- function(y, x, lambda,
                         grid = seq(0, sd(y, na.rm = TRUE)/mean(diff(x))*10, sd(y, na.rm = TRUE)/mean(diff(x))/10),
                         ic = c("bic", "ebic", "hq", "aic", "aicc"),
                         max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                         ebic_xi = 1
                         ) {
  ic <- match.arg(ic)
  last_ic     <- Inf
  n_increases <- 0L
  gamma_init  <- numeric(0)
  for (M in grid) {
    out <- ssj_atw(y = y, x = x, lambda = lambda,
                   maxsum = M, max_iter = max_iter, tol = tol,
                   learning_rate = learning_rate, gamma_init = gamma_init,
                   ebic_xi = ebic_xi)
    gamma_init <- out$gamma   # warm start for the next grid point
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
