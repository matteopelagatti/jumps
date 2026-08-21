#' Smoothing splines with automatic jump detection using the alternating
#' time-warping algorithm
#'
#' Same model as \code{\link{ssj_fix}}, but fitted with the alternating
#' time-warping (ATW) algorithm: a two-step block coordinate descent in which
#' the \emph{f}-step solves a banded linear system in \eqn{O(n)} operations and
#' the \eqn{\gamma}-step performs a gradient descent on the roughness energy
#' followed by a projection onto the simplex \eqn{\{\gamma \ge 0, \sum_i
#' \gamma_i \le M\}}.
#'
#' Internally \code{x} is mapped to \eqn{[0,1]} and \code{y} is scaled to unit
#' variance, so that \code{lambda} and \code{maxsum} are interpreted in the
#' original units of the data. Tied values of \code{x} are collapsed by
#' averaging the corresponding \code{y}; the fitted values are therefore
#' returned at the (possibly shorter) vector of unique, sorted \code{x} values
#' available in the \code{x} slot of the output.
#'
#' @param y numeric vector with the response (no missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param lambda numeric scalar with the smoothing parameter.
#' @param maxsum maximum sum of the \eqn{\gamma_t}, that is the warping budget
#'   \eqn{M} of the constraint \eqn{\sum_t \gamma_t \le M}. It is expressed in
#'   the units of \code{x} and defaults to the range of \code{x}, the same
#'   upper end as the grid searched by \code{\link{auto_ssj_atw}}.
#' @param max_iter maximum number of iterations (default 100).
#' @param tol tolerance for numerical convergence (default 1e-4).
#' @param learning_rate initial and fallback step size for the Barzilai-Borwein
#'   adaptive scheme (default 0.1).
#' @param gamma_init optional numeric vector of length \code{n-1} with a
#'   warm start for \eqn{\gamma}; ignored when its length does not match
#'   \code{n-1}.
#' @param weights optional numeric vector of positive observation weights, one
#'   per element of \code{y}. The data term becomes
#'   \eqn{\sum_i w_i (y_i - f_i)^2}, so a point with weight \eqn{w_i} counts as
#'   the mean of \eqn{w_i} independent observations, and the information
#'   criteria use \eqn{\sum_i w_i} in place of the number of design points.
#'   Weights are averaged along with \code{y} when tied values of \code{x} are
#'   collapsed. The default, \code{NULL}, means unit weights.
#' @param ebic_xi numeric scalar in \code{[0,1]}: strength of the extended-BIC
#'   penalty for searching over the \code{n-1} candidate jump locations
#'   (default 0.25; 0 reduces \code{ebic} to plain \code{bic}). The default is
#'   calibrated rather than taken from theory; see the note below.
#'
#' @section Choosing \code{ebic_xi}:
#' The default is the value that minimises worst-case regret in root mean
#' squared error over the 48-cell design of
#' \code{inst/scripts/xi_calibration.R}. Unlike \code{\link{ssj_mle}} and
#' \code{\link{ssj_kfs}}, whose optima are flat in \eqn{n}, the best value for
#' this estimator falls steadily with the sample size -- about 0.32, 0.22, 0.15
#' and 0.12 at \eqn{n} = 50, 100, 200 and 400. What is constant instead is the
#' total charge per declared discontinuity,
#' \eqn{s\log n + 2\xi\log(n-1)} with \eqn{s \approx 0.42} the degrees of
#' freedom each activated warp adds, which held at 3.9 across the design
#' (trend \eqn{-0.10} per \eqn{\log n}, 95\% interval \eqn{[-0.23, 0.02]}).
#' That gives the rule
#' \deqn{\xi(n) = \frac{3.93 - 0.419\log n}{2\log(n-1)},}
#' which is worth using in place of the constant when \eqn{n} is far from a
#' few hundred: at \eqn{n = 400} it improved root mean squared error by 4.4\%
#' over the fixed 0.25. The constant is the default because it is the value
#' that behaves acceptably everywhere, and because a default that changes with
#' the data is easy to overlook.
#'
#' Note that the charge of about 3.9 on the \eqn{-2\log L} scale is close to 4,
#' that is roughly two log-likelihood units per discontinuity -- what AIC
#' charges for one free parameter.
#'
#' @returns An object of class \code{ssj} (a list, so that all the slots below
#'   remain directly accessible) with a \code{\link[=print.ssj]{print}} and a
#'   \code{\link[=plot.ssj]{plot}} method. Its slots are:
#' \itemize{
#'  \item \code{nobs}: number of distinct sampling sites used.
#'  \item \code{df}: model's effective degrees of freedom.
#'  \item \code{maxsum}: the \code{maxsum} argument.
#'  \item \code{pars}: named vector (sigma, lambda), where \code{sigma} is the
#'    unbiased residual standard deviation \eqn{\sqrt{rss/(n - edf)}} that the
#'    information criteria are evaluated at, and \code{lambda} is the argument,
#'    recorded so that it can be recovered from a fit.
#'  \item \code{loglik}: Gaussian log-likelihood evaluated at the unbiased
#'    residual variance \eqn{\sigma^2 = rss/(n - edf)}, on the scale of
#'    \code{y} as supplied.
#'  \item \code{gcv}: generalized cross-validation criterion.
#'  \item \code{gamma}: vector of additional times (time warps) of length
#'    \code{nobs - 1}.
#'  \item \code{ic}: vector of information criteria (aic, aicc, bic, hq, ebic).
#'  \item \code{smoothed_level}: vector with the smoothing spline with jumps.
#'  \item \code{x}: the ordered vector of distinct \code{x} values.
#' }
#'
#' @references
#' Chen, J. and Chen, Z. (2008). Extended Bayesian information criteria for
#' model selection with large model spaces. \emph{Biometrika}, \strong{95}(3),
#' 759--771. \doi{10.1093/biomet/asn034}
#'
#' @seealso \code{\link{auto_ssj_atw}}, \code{\link{ssj_fix}},
#'   \code{\link{ssj_kfs}}, \code{\link{ssj_mle}}.
#'
#' @examples
#' # 'maxsum' is a warping budget in the units of x; the range of eruptions
#' # is about 3.5, so a budget of 1 is already generous.
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- ssj_atw(y = faithful$waiting,
#'                   x = faithful$eruptions,
#'                   lambda = 1,
#'                   maxsum = 1)
#' lines(x = of_jum$x, y = of_jum$smoothed_level, col = "red", lwd = 2)
#'
#' @export
ssj_atw <- function(y, x, lambda, maxsum = max(x) - min(x), weights = NULL,
                    max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                    gamma_init = numeric(0), ebic_xi = 0.25) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  wt <- ssj_check_weights(weights, length(y))
  if (length(y) != length(x)) stop("'x' and 'y' must have the same length")
  if (anyNA(x)) stop("'x' must not contain missing values")
  ord <- order(x)
  y <- y[ord]
  x <- x[ord]
  wt <- wt[ord]
  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    # collapse ties: the mean of the responses, carrying the summed weight, is
    # the sufficient statistic for the group
    y   <- as.numeric(tapply(y * wt, grp, sum) / tapply(wt, grp, sum))
    wt  <- as.numeric(tapply(wt, grp, sum))
    x   <- x[!duplicated(x)]
  }
  # The banded solver of the f-step has no missing-value handling: NAs would
  # propagate to the whole fit.  Use ssj_kfs() or ssj_mle() when y has gaps.
  if (anyNA(y)) {
    stop("'y' must not contain missing values: use ssj_kfs() or ssj_mle() instead")
  }
  n       <- length(x)
  # Force the default of maxsum while x is still on its original scale (here
  # the normalised copy is a new object, x_w, so this is belt and braces).
  force(maxsum)
  # Normalise: x -> [0,1], y -> unit variance.
  # Cubic-spline penalty lambda * integral(f'')^2 dx => lambda_norm = lambda / x_range^3.
  x_range  <- x[n] - x[1L]
  x_min_s  <- x[1L]
  y_scale  <- sqrt(var(y))
  if (!is.finite(y_scale) || y_scale < .Machine$double.eps) y_scale <- 1
  if (x_range < .Machine$double.eps) stop("'x' has zero range")
  maxsum_orig  <- maxsum
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
                                ebic_xi = ebic_xi, w_in = wt)
  # Back to the scale of the data as supplied.  The C++ routine sees y already
  # divided by y_scale, so its residual variance is smaller by y_scale^2 and
  # its log-likelihood is larger by wsum*log(y_scale); the criteria, which
  # enter as -2*loglik, are correspondingly smaller by 2*wsum*log(y_scale).
  # (gcv is already corrected above by the y_scale^2 factor.)  The offset is
  # the same at every value of maxsum, so it cancels in the comparisons that
  # auto_ssj_atw() makes and cannot change which model is selected.
  wsum   <- sum(wt)
  ll_adj <- wsum * log(y_scale)
  # sigma is the unbiased residual standard deviation the criteria are built
  # on, back on the scale of y; lambda is reported so that the smoothing
  # constant can be recovered from a fit, as it can from the other estimators.
  sigma_hat <- if (n - out$edf > 0) sqrt(out$rss / (n - out$edf)) * y_scale
               else NA_real_
  res <- list(nobs          = n,
       df            = out$edf,
       maxsum        = maxsum_orig,
       pars          = c(sigma = sigma_hat, lambda = lambda),
       loglik        = out$loglik - ll_adj,
       gcv           = out$gcv * y_scale^2,
       gamma         = out$gamma * x_range,
       ic            = c(aic  = out$aic  + 2 * ll_adj,
                         aicc = out$aicc + 2 * ll_adj,
                         bic  = out$bic  + 2 * ll_adj,
                         hq   = out$hq   + 2 * ll_adj,
                         ebic = out$ebic + 2 * ll_adj),
       smoothed_level = out$f * y_scale,
       smoothed_slope = out$df * y_scale / x_range,
       x              = x
  )
  new_ssj(res, method = "atw", y = y, call = match.call())
}

#' Automatic selection of the optimal smoothing spline with jumps using the
#' alternating time-warping algorithm
#'
#' The regularization parameter of the ATW algorithm is the maximal sum of
#' additional times. This value has to be passed to the \code{\link{ssj_atw}}
#' function. \code{auto_ssj_atw} runs \code{ssj_atw} on a grid of regularization
#' constants and returns the output of \code{ssj_atw} selected by the chosen
#' information criterion. Each grid point is warm-started from the solution at
#' the previous point, and the search stops early when the selected criterion
#' increases for three consecutive grid values.
#'
#' @param y numeric vector with the response (no missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param lambda numeric scalar with the smoothing parameter.
#' @param grid numeric vector of \code{maxsum} values (in the units of \code{x})
#'   to search; the default is 21 equally spaced values from 0 to the range of
#'   \code{x}.
#' @param ic string with the information criterion used for the choice: the
#'   default is \code{"ebic"} (recommended when jumps are expected to be rare
#'   relative to the \code{n-1} candidate locations, see \code{ebic_xi}), but
#'   \code{"bic"}, \code{"hq"}, \code{"aic"} and \code{"aicc"} are also
#'   available.
#' @param max_iter maximum number of iterations of \code{\link{ssj_atw}}
#'   (default 100).
#' @param tol tolerance for the convergence of \code{\link{ssj_atw}}
#'   (default 1e-4).
#' @param learning_rate learning rate of \code{\link{ssj_atw}} (default 0.1).
#' @param weights optional numeric vector of positive observation weights,
#'   passed to \code{\link{ssj_atw}}; see its help page.
#' @param ebic_xi numeric scalar in \code{[0,1]} passed to
#'   \code{\link{ssj_atw}}: strength of the extended-BIC search penalty
#'   (default 0.25; only relevant when \code{ic = "ebic"}). See the
#'   "Choosing ebic_xi" section of \code{\link{ssj_atw}} for how it was
#'   calibrated, why it differs from the likelihood-based estimators, and for
#'   the rule that adapts it to the sample size.
#'
#' @returns The output of the \code{\link{ssj_atw}} function corresponding to
#'   the best choice according to the selected information criterion.
#'
#' @seealso \code{\link{ssj_atw}}.
#'
#' @examples
#' plot(faithful$eruptions, faithful$waiting)
#' of_jum <- auto_ssj_atw(y = faithful$waiting,
#'                        x = faithful$eruptions,
#'                        lambda = 1)
#' lines(x = of_jum$x, y = of_jum$smoothed_level, col = "red", lwd = 2)
#'
#' @export
auto_ssj_atw <- function(y, x, lambda, weights = NULL,
                         grid = seq(0, max(x) - min(x), length.out = 21L),
                         ic = c("ebic", "bic", "hq", "aic", "aicc"),
                         max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                         ebic_xi = 0.25
                         ) {
  ic <- match.arg(ic)
  last_ic     <- Inf
  n_increases <- 0L
  gamma_init  <- numeric(0)
  for (M in grid) {
    out <- ssj_atw(y = y, x = x, lambda = lambda, weights = weights,
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
  best$call <- match.call()
  best
}
