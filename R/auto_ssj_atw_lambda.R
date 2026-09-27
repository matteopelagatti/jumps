# Selection of the smoothing constant for the alternating time-warping
# algorithm, without recourse to ssj_mle().
#
# Both routines split the problem in two: choose lambda first, then hand it to
# auto_ssj_atw(), which selects the warping budget M by information criterion
# exactly as before.  Choosing the two jointly is deliberately avoided; see the
# note in auto_ssj_atw_gcv() on why a criterion minimised over lambda AND M
# collapses onto an interpolating fit.

# Default grid of smoothing constants.
#
# lambda multiplies the integrated squared second derivative, so it has the
# units of x^3 and a grid fixed in absolute terms would only suit one scaling
# of x.  ssj_atw() normalises x to [0, 1] internally, and it is that normalised
# lambda which is scale free, so the grid is built there and mapped back.
#
# The lower end is 1e-6 rather than 0.  Generalized cross-validation is not
# reliable at arbitrarily small lambda here: the criterion keeps improving as
# the fit approaches interpolation, and without a floor it selects the smallest
# value on offer.  1e-6 is where that behaviour set in on the calibration
# design; the selected value is checked against both ends of the grid and a
# warning is issued if it lands there.
ssj_lambda_grid <- function(x, n_lambda = 9L, from = 1e-6, to = 1e-2) {
  x_range <- max(x) - min(x)
  if (!is.finite(x_range) || x_range <= 0) stop("'x' has zero range")
  10^seq(log10(from), log10(to), length.out = n_lambda) * x_range^3
}

ssj_warn_boundary <- function(lambda_grid, chosen, what) {
  if (length(lambda_grid) < 2L) return(invisible(NULL))
  if (isTRUE(all.equal(chosen, min(lambda_grid))))
    warning(what, " selected the smallest lambda on the grid (", signif(chosen, 3),
            "): the criterion is still improving there, so the value is a ",
            "boundary solution rather than a minimum. Widen 'lambda_grid' ",
            "downwards only if you are prepared for a near-interpolating fit.",
            call. = FALSE)
  else if (isTRUE(all.equal(chosen, max(lambda_grid))))
    warning(what, " selected the largest lambda on the grid (", signif(chosen, 3),
            "): extend 'lambda_grid' upwards.", call. = FALSE)
  invisible(NULL)
}

#' Smoothing splines with jumps by alternating time warping, with the smoothing
#' constant chosen by generalized cross-validation
#'
#' Selects \code{lambda} by minimising the generalized cross-validation
#' criterion of a plain smoothing spline -- that is, of \code{\link{ssj_atw}}
#' with the warping budget held at zero -- and then passes it to
#' \code{\link{auto_ssj_atw}}, which chooses the budget \eqn{M} by information
#' criterion in the usual way. Unlike the two-step recipe in which
#' \code{lambda} comes from \code{\link{ssj_mle}}, nothing here depends on the
#' likelihood-based estimator, so the alternating algorithm stands on its own.
#'
#' The cost is one fit per candidate \code{lambda}, because \code{ssj_atw}
#' already returns the exact effective degrees of freedom, the trace of the
#' smoothing matrix, so the criterion is available in closed form and no
#' resampling is needed. That makes this route roughly two orders of magnitude
#' cheaper than \code{\link{auto_ssj_atw_cv}}.
#'
#' \strong{Why the budget is held at zero.} Minimising the criterion over
#' \code{lambda} and \eqn{M} jointly does not work. The warps give the fit a
#' way to interpolate: the effective degrees of freedom approach the number of
#' observations, the denominator of the criterion collapses, and the smallest
#' \code{lambda} on offer wins. On the calibration design that cost 58\% in root
#' mean squared error and reduced the Matthews correlation of the detected
#' discontinuities from 0.54 to 0.08. With no budget there are no warps to
#' exploit and the classical behaviour is restored.
#'
#' \strong{Know what this routine can and cannot do.} At zero budget the spline
#' has no way to represent a discontinuity except by bending sharply, so when
#' the data contain clear jumps and little noise the criterion keeps improving
#' as \code{lambda} falls and the selection lands on the floor of the grid.
#' On a design with three jumps this happened in 83\% of samples at
#' \eqn{n = 100} and in every sample at \eqn{n = 200} when \eqn{\sigma = 0.05},
#' against 7\% and 23\% at \eqn{\sigma = 0.15} and about 1\% with no jumps at
#' all. In that regime the routine is not choosing \code{lambda}; it is
#' returning the floor, and the fit is only as good as that floor happens to
#' be. The resulting accuracy is in fact competitive -- root mean squared error
#' within a few percent of the alternatives -- but the value is not a minimum
#' of the criterion, and the warning issued in that case should be taken at
#' face value.
#'
#' \code{\link{auto_ssj_atw_cv}} does not have this failure mode: its floor
#' bound in 0 to 23\% of the same samples and the \code{lambda} it chooses
#' tracks the maximum-likelihood value closely. Prefer it when the cost is
#' affordable, and prefer this routine when it is not, or when the signal is
#' noisy enough that the degeneracy does not bite.
#'
#' @param y numeric vector with the response (no missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param weights optional numeric vector of positive observation weights,
#'   passed to \code{\link{ssj_atw}}; see its help page.
#' @param lambda_grid numeric vector of candidate smoothing constants, in the
#'   units of \code{x}. The default is nine values spanning
#'   \code{1e-6} to \code{1e-2} once \code{x} has been scaled to unit range,
#'   which is the range over which the criterion behaves; values below the
#'   floor make it degenerate. A warning is issued when the selected value is
#'   an endpoint of the grid.
#' @param grid numeric vector of \code{maxsum} values passed to
#'   \code{\link{auto_ssj_atw}}; the default is 21 equally spaced values from
#'   0 to the range of \code{x}.
#' @param ic string with the information criterion used to choose \eqn{M},
#'   passed to \code{\link{auto_ssj_atw}} (default \code{"ebic"}).
#' @param max_iter,tol,learning_rate,ebic_xi passed to \code{\link{ssj_atw}}.
#'
#' @returns The object returned by \code{\link{auto_ssj_atw}} at the selected
#'   \code{lambda}, of class \code{ssj}, with two extra slots:
#' \itemize{
#'   \item \code{lambda_path}: data frame of the candidate \code{lambda} values
#'     and their generalized cross-validation criterion.
#'   \item \code{lambda_source}: the string \code{"gcv"}.
#' }
#'
#' @references
#' Craven, P. and Wahba, G. (1978). Smoothing noisy data with spline functions.
#' \emph{Numerische Mathematik}, \strong{31}(4), 377--403.
#' \doi{10.1007/BF01404567}
#'
#' @seealso \code{\link{auto_ssj_atw_cv}}, \code{\link{auto_ssj_atw}},
#'   \code{\link{ssj_atw}}.
#'
#' @examples
#' set.seed(11)
#' x <- sort(runif(80))
#' y <- sin(2*pi*x) + 0.6*(x > 0.45) + rnorm(80, sd = 0.1)
#' fit <- auto_ssj_atw_gcv(y, x)
#' plot(x, y, pch = 16, cex = 0.6, col = "grey60")
#' lines(fit$x, fit$smoothed_level, col = "firebrick", lwd = 2)
#' fit$pars["lambda"]
#'
#' @export
auto_ssj_atw_gcv <- function(y, x, weights = NULL, lambda_grid = NULL,
                             grid = seq(0, max(x) - min(x), length.out = 21L),
                             ic = c("ebic", "bic", "hq", "aic", "aicc"),
                             max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                             ebic_xi = 0.25) {
  ic <- match.arg(ic)
  y  <- as.numeric(y); x <- as.numeric(x)
  if (length(y) != length(x)) stop("'x' and 'y' must have the same length")
  if (anyNA(x) || anyNA(y))
    stop("'x' and 'y' must not contain missing values: ssj_atw() has no ",
         "missing-value handling")
  if (is.null(lambda_grid)) lambda_grid <- ssj_lambda_grid(x)
  lambda_grid <- sort(as.numeric(lambda_grid))
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0))
    stop("'lambda_grid' must be a vector of positive, finite values")

  # one fit per lambda, at zero budget
  gcv <- vapply(lambda_grid, function(l) {
    f <- tryCatch(ssj_atw(y = y, x = x, lambda = l, maxsum = 0,
                          weights = weights, max_iter = max_iter, tol = tol,
                          learning_rate = learning_rate, ebic_xi = ebic_xi),
                  error = function(e) NULL)
    if (is.null(f) || !is.finite(f$gcv)) Inf else f$gcv
  }, numeric(1))
  if (all(!is.finite(gcv)))
    stop("every candidate lambda failed: check 'lambda_grid' and the data")

  lambda_hat <- lambda_grid[which.min(gcv)]
  ssj_warn_boundary(lambda_grid, lambda_hat, "auto_ssj_atw_gcv()")

  out <- auto_ssj_atw(y = y, x = x, lambda = lambda_hat, weights = weights,
                      grid = grid, ic = ic, max_iter = max_iter, tol = tol,
                      learning_rate = learning_rate, ebic_xi = ebic_xi)
  out$lambda_path   <- data.frame(lambda = lambda_grid, gcv = gcv)
  out$lambda_source <- "gcv"
  out$call          <- match.call()
  out
}

#' Smoothing splines with jumps by alternating time warping, with the smoothing
#' constant chosen by cross-validation
#'
#' Selects \code{lambda} by \eqn{K}-fold cross-validation over the smoothing
#' constant and the warping budget jointly, keeps \emph{only} the selected
#' \code{lambda}, and passes it to \code{\link{auto_ssj_atw}}, which chooses
#' the budget \eqn{M} by information criterion. Like
#' \code{\link{auto_ssj_atw_gcv}} it makes the alternating algorithm
#' self-contained; it is markedly more expensive but does not need the floor on
#' the grid that generalized cross-validation requires.
#'
#' \strong{Why the budget is re-selected afterwards.} The cross-validated
#' budget is not used. Cross-validation estimates prediction error, and
#' prediction error is nearly indifferent to a handful of spurious warps: on
#' the calibration design, selecting \eqn{M} this way declared substantially
#' more discontinuities than the information criterion at the same accuracy.
#' The smoothing constant is a different matter -- it governs a global
#' bias-variance trade-off to which every held-out point contributes -- and
#' there cross-validation is competitive with maximum likelihood. Hence the
#' split: cross-validate \code{lambda}, leave \eqn{M} to the criterion.
#'
#' Tied values of \code{x} require a choice about the partition. Folds are drawn
#' over the distinct values of \code{x}, so every replicate at a tied covariate
#' value is held out together. Drawing them over the observations instead would
#' leave replicates of a held-out site in the training sample, where
#' \code{\link{ssj_atw}} collapses them and pins the fit at that very covariate
#' value; the criterion would then estimate the error of predicting a fresh
#' observation at an already observed \code{x} rather than at a new one, and is
#' optimistic for the present purpose -- by roughly a third under heavy
#' replication. For a design with no ties the two partitions coincide.
#'
#' @param y numeric vector with the response (no missing values allowed).
#' @param x numeric vector with the covariate (no missing values allowed).
#' @param weights optional numeric vector of positive observation weights,
#'   passed to \code{\link{ssj_atw}} and used to weight the held-out squared
#'   errors; see its help page.
#' @param lambda_grid numeric vector of candidate smoothing constants, in the
#'   units of \code{x}; the default is as in \code{\link{auto_ssj_atw_gcv}}.
#'   No floor is needed here, so it may be extended downwards freely.
#' @param nfolds integer number of cross-validation folds (default 5). Folds
#'   are drawn over the \emph{distinct} values of \code{x}, so all replicates
#'   at a tied covariate value are held out together; \code{nfolds} must
#'   therefore not exceed the number of distinct values. See Details.
#' @param folds optional integer vector of fold labels, one per observation,
#'   for reproducibility or for grouped designs; overrides \code{nfolds}. A
#'   partition that puts replicates of the same \code{x} in different folds
#'   draws a warning, since the resulting criterion is optimistic.
#' @param grid numeric vector of \code{maxsum} values, used both for the
#'   cross-validation surface and for the final selection of \eqn{M}.
#' @param ic string with the information criterion used to choose \eqn{M},
#'   passed to \code{\link{auto_ssj_atw}} (default \code{"ebic"}).
#' @param max_iter,tol,learning_rate,ebic_xi passed to \code{\link{ssj_atw}}.
#'
#' @returns The object returned by \code{\link{auto_ssj_atw}} at the selected
#'   \code{lambda}, of class \code{ssj}, with three extra slots:
#' \itemize{
#'   \item \code{lambda_path}: data frame of the candidate \code{lambda} values
#'     and their best cross-validated mean squared prediction error.
#'   \item \code{cv_surface}: the full \code{lambda} by \code{maxsum} matrix of
#'     cross-validated errors.
#'   \item \code{lambda_source}: the string \code{"cv"}.
#' }
#'
#' @seealso \code{\link{auto_ssj_atw_gcv}}, \code{\link{auto_ssj_atw}},
#'   \code{\link{ssj_atw}}.
#'
#' @examples
#' set.seed(11)
#' x <- sort(runif(60))
#' y <- sin(2*pi*x) + 0.6*(x > 0.45) + rnorm(60, sd = 0.1)
#' # a short grid keeps the example quick
#' fit <- auto_ssj_atw_cv(y, x, nfolds = 3,
#'                        lambda_grid = 10^seq(-6, -3, length.out = 4),
#'                        grid = seq(0, 1, length.out = 6))
#' plot(x, y, pch = 16, cex = 0.6, col = "grey60")
#' lines(fit$x, fit$smoothed_level, col = "steelblue", lwd = 2)
#'
#' @export
auto_ssj_atw_cv <- function(y, x, weights = NULL, lambda_grid = NULL,
                            nfolds = 5L, folds = NULL,
                            grid = seq(0, max(x) - min(x), length.out = 21L),
                            ic = c("ebic", "bic", "hq", "aic", "aicc"),
                            max_iter = 100, tol = 1e-4, learning_rate = 0.1,
                            ebic_xi = 0.25) {
  ic <- match.arg(ic)
  y  <- as.numeric(y); x <- as.numeric(x)
  if (length(y) != length(x)) stop("'x' and 'y' must have the same length")
  if (anyNA(x) || anyNA(y))
    stop("'x' and 'y' must not contain missing values: ssj_atw() has no ",
         "missing-value handling")
  n  <- length(y)
  wt <- ssj_check_weights(weights, n)
  if (is.null(lambda_grid)) lambda_grid <- ssj_lambda_grid(x)
  lambda_grid <- sort(as.numeric(lambda_grid))
  if (!length(lambda_grid) || any(!is.finite(lambda_grid)) || any(lambda_grid <= 0))
    stop("'lambda_grid' must be a vector of positive, finite values")
  grid <- sort(as.numeric(grid))

  # Folds are drawn over the DISTINCT values of x, not over the observations.
  # ssj_atw() collapses tied x in whatever sample it is given, so a partition
  # made observation by observation leaves replicates of a held-out site in the
  # training set, and the held-out point is then scored at a covariate value the
  # training sample already pins down.  The criterion that results estimates the
  # error of predicting a fresh observation at an ALREADY OBSERVED x, which is
  # optimistic for the quantity wanted here -- the error at a new x -- by about
  # a third under heavy replication.  Grouping by site targets the right
  # quantity.  With no ties the two partitions coincide, so nothing changes for
  # a design with distinct covariate values.
  site   <- match(x, sort(unique(x)))
  nsite  <- max(site)
  if (is.null(folds)) {
    nfolds <- as.integer(nfolds)
    if (!is.finite(nfolds) || nfolds < 2L || nfolds > nsite)
      stop("'nfolds' must be between 2 and the number of distinct values of ",
           "'x' (", nsite, ")")
    folds <- sample(rep_len(seq_len(nfolds), nsite))[site]
  } else {
    folds <- as.integer(folds)
    if (length(folds) != n) stop("'folds' must have the same length as 'y'")
    if (length(unique(folds)) < 2L) stop("'folds' must define at least 2 folds")
    # A supplied partition is respected -- the caller may have a grouping of
    # their own -- but a split site is almost always an oversight, so say so.
    split_sites <- sum(tapply(folds, site, function(z) length(unique(z))) > 1L)
    if (split_sites > 0L)
      warning("'folds' splits the replicates of ", split_sites,
              " tied value(s) of 'x' across folds: the cross-validated error ",
              "is then optimistic, because a held-out point is scored at a ",
              "covariate value that is also in the training sample.",
              call. = FALSE)
  }

  err <- den <- matrix(0, length(lambda_grid), length(grid),
                       dimnames = list(signif(lambda_grid, 4), signif(grid, 4)))
  for (k in sort(unique(folds))) {
    tr <- folds != k; te <- !tr
    if (!any(te) || sum(tr) < 4L) next
    for (a in seq_along(lambda_grid)) {
      gi <- numeric(0)                       # warm start along the M grid
      for (j in seq_along(grid)) {
        f <- tryCatch(ssj_atw(y = y[tr], x = x[tr], lambda = lambda_grid[a],
                              maxsum = grid[j], weights = wt[tr],
                              max_iter = max_iter, tol = tol,
                              learning_rate = learning_rate,
                              gamma_init = gi, ebic_xi = ebic_xi),
                      error = function(e) NULL)
        if (is.null(f)) { gi <- numeric(0); next }
        gi <- f$gamma
        # held-out points outside the training range are clamped rather than
        # extrapolated: beyond the data the fit is only linearly extended, and
        # a couple of edge points would otherwise dominate the fold error
        xe <- pmin(pmax(x[te], min(f$x)), max(f$x))
        pr <- tryCatch(stats::predict(f, xe), error = function(e) NULL)
        if (is.null(pr) || anyNA(pr)) next
        err[a, j] <- err[a, j] + sum(wt[te] * (y[te] - pr)^2)
        den[a, j] <- den[a, j] + sum(wt[te])
      }
    }
  }
  cv <- ifelse(den > 0, err/den, Inf)
  if (all(!is.finite(cv)))
    stop("every candidate (lambda, maxsum) pair failed during cross-validation")

  # keep the lambda of the joint minimiser; the budget is re-selected below
  lambda_hat <- lambda_grid[which(cv == min(cv), arr.ind = TRUE)[1L, "row"]]
  ssj_warn_boundary(lambda_grid, lambda_hat, "auto_ssj_atw_cv()")

  out <- auto_ssj_atw(y = y, x = x, lambda = lambda_hat, weights = weights,
                      grid = grid, ic = ic, max_iter = max_iter, tol = tol,
                      learning_rate = learning_rate, ebic_xi = ebic_xi)
  out$lambda_path   <- data.frame(lambda = lambda_grid,
                                  cv = apply(cv, 1L, min))
  out$cv_surface    <- cv
  out$lambda_source <- "cv"
  out$call          <- match.call()
  out
}
