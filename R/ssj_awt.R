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
ssj_atw <- function(y, x, lambda, maxsum = max(x) - min(x),
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
  n       <- length(x)
  # Normalise: x -> [0,1], y -> unit variance.
  # Cubic-spline penalty lambda * integral(f'')^2 dx => lambda_norm = lambda / x_range^3.
  x_range  <- x[n] - x[1L]
  x_min_s  <- x[1L]
  y_scale  <- sqrt(var(y, na.rm = TRUE))
  if (!is.finite(y_scale) || y_scale < .Machine$double.eps) y_scale <- 1
  if (x_range < .Machine$double.eps) stop("'x' has zero range")
  maxsum_orig  <- maxsum
  lambda_orig  <- lambda
  x_w          <- (x - x_min_s) / x_range
  y_w          <- y / y_scale
  lambda_n     <- lambda / x_range^3
  maxsum_n     <- maxsum / x_range
  gamma_init_n <- if (length(gamma_init) == n - 1L) gamma_init / x_range else numeric(0)
  out <- solve_jump_spline_fast(x_in = x_w, y_in = y_w,
                                lambda = lambda_n, M = maxsum_n,
                                max_iter = max_iter, tol = tol,
                                learning_rate = learning_rate,
                                gamma_init = gamma_init_n,
                                ebic_xi = ebic_xi)
  list(nobs          = n,
       df            = out$edf,
       maxsum        = maxsum_orig,
       loglik        = out$loglik,
       gcv           = out$gcv * y_scale^2,
       gamma         = out$gamma * x_range,
       ic            = c(aic  = out$aic,
                         aicc = out$aicc,
                         bic  = out$bic,
                         hq   = out$hq,
                         ebic = out$ebic),
       smoothed_level = out$f * y_scale,
       x              = x
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
#' @param grid numeric vector of maxsum values (in x units) to search; default
#' is 21 equally spaced values from 0 to the x range;
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
                         grid = seq(0, max(x) - min(x), length.out = 21L),
                         ic = c("ebic", "bic", "hq", "aic", "aicc"),
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

#' Joint selection of smoothing parameter and jump budget for ATW splines
#'
#' Extends \code{auto_ssj_atw} by also searching over a grid of smoothing
#' parameters \code{lambda}. For each candidate \code{lambda}, the optimal
#' jump budget (\code{maxsum}) is chosen by the information criterion \code{ic}
#' via \code{auto_ssj_atw}; the best \code{lambda} is then selected by
#' generalised cross-validation (GCV). This two-criterion strategy is
#' theoretically motivated: GCV targets the continuous smoothing dimension,
#' while the discrete jump-location selection is handled by a sparse IC.
#'
#' @param y vector with the y (missing values allowed).
#' @param x vector with the x (no missing values allowed).
#' @param lambda_grid numeric vector of candidate smoothing parameters (in the
#'   same units as the user-supplied lambda, i.e.\ original \eqn{x} scale).
#'   Defaults to 20 log-spaced values spanning the normalised-coordinate GCV
#'   window \eqn{[10^0, 10^{10}]} scaled back by \eqn{x\_range^3} so the
#'   search is always in the right range regardless of \eqn{x} units.
#' @param maxsum_grid numeric vector of candidate jump budgets (in x units)
#'   passed to \code{auto_ssj_atw} for each \code{lambda} evaluation; default
#'   is 21 equally spaced values from 0 to the x range.
#' @param ic string with the information criterion used to select
#'   \code{maxsum} within each \code{lambda}: "ebic" (default, recommended for
#'   sparse jumps), "bic", "hq", "aic", or "aicc".
#' @param max_iter maximum number of iterations in \code{ssj_atw} (default 100).
#' @param tol tolerance for convergence of \code{ssj_atw} (default 1e-4).
#' @param learning_rate learning rate of \code{ssj_atw} (default 0.1).
#' @param ebic_xi numeric scalar in [0,1] passed to \code{ssj_atw}: strength
#'   of the Extended-BIC search penalty (default 1; only relevant when
#'   \code{ic = "ebic"}).
#'
#' @returns The output of \code{ssj_atw} for the selected
#'   \code{(lambda, maxsum)} pair, with an additional field \code{lambda}
#'   recording the chosen smoothing parameter.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- auto_ssj_atw_gcv(y = faithful$waiting,
#'                             x = faithful$eruptions)
#' lines(x = sort(faithful$eruptions),
#'       y = of_jum$smoothed_level, col = "red", lwd = 2)
#'
#' @export
auto_ssj_atw_gcv <- function(y, x,
                              lambda_grid = NULL,
                              maxsum_grid = seq(0, max(x) - min(x), length.out = 21L),
                              ic = c("ebic", "bic", "hq", "aic", "aicc"),
                              max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                              ebic_xi = 1) {
  ic <- match.arg(ic)

  if (is.null(lambda_grid)) {
    # ssj_atw normalises x -> [0,1] internally, scaling lambda by 1/x_range^3.
    # For unit-variance signals on [0,1] the GCV optimum is roughly in
    # [10^2, 10^8]; supplying lambda in original units requires multiplying by
    # x_range^3 so the normalised values stay in that target window.
    x_range_local <- max(x) - min(x)
    lambda_grid   <- 10^seq(0, 10, length.out = 20) * x_range_local^3
  }

  best_gcv    <- Inf
  best_out    <- NULL
  best_lambda <- NA_real_

  for (lam in lambda_grid) {
    out <- auto_ssj_atw(y = y, x = x, lambda = lam,
                        grid = maxsum_grid, ic = ic,
                        max_iter = max_iter, tol = tol,
                        learning_rate = learning_rate,
                        ebic_xi = ebic_xi)
    if (out$gcv < best_gcv) {
      best_gcv    <- out$gcv
      best_out    <- out
      best_lambda <- lam
    }
  }

  best_out$lambda <- best_lambda
  best_out
}
