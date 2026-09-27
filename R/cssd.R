#' Total-variation denoising for irregularly sampled signals
#'
#' Solves the total-variation denoising problem
#' \deqn{\min_u \frac{1}{2}\|y - u\|_2^2 + \lambda\|Du\|_1,}
#' where \eqn{D} is the unweighted first-difference operator
#' \eqn{(Du)_i = u_{i+1} - u_i}, using ADMM with a pre-factorized
#' \eqn{I + \rho D^\top D} linear system. The minimizer \eqn{u} is piecewise
#' constant, with non-zero first differences at the detected jump locations.
#'
#' This implements the simplified total-variation (first-order) variant related
#' to the CSSD framework of Storath and Weinmann (2024). The full CSSD
#' (piecewise cubic smoothing spline with discontinuities) is implemented in
#' \code{\link{ssj_cssd_dp}}; the total-variation formulation used here
#' replaces the \eqn{C^2} roughness penalty with an \eqn{\ell_1} penalty on
#' first differences and is therefore much faster, at the cost of a piecewise
#' constant (rather than piecewise cubic) fit.
#'
#' The unweighted \eqn{D} (entries \eqn{\pm 1}) is used rather than the
#' slope-based \eqn{(Du)_i = (u_{i+1}-u_i)/\delta_i}: the \eqn{1/\delta}
#' weighting inflates the eigenvalues of \eqn{D^\top D} by \eqn{1/h^2 \sim
#' 10^4} for \eqn{h\sim 0.01}, which drives the ADMM \eqn{u}-update towards
#' zero in the first iteration and stalls it permanently.
#'
#' @param y numeric vector of observations (no missing values).
#' @param x numeric vector of sampling locations (no missing values; it need
#'   not be sorted). If \code{NULL} (default), regular spacing \code{1:n} is
#'   assumed.
#' @param lambda positive regularization parameter controlling the
#'   soft-threshold \eqn{\lambda/\rho} applied to the first differences of
#'   \eqn{u}. Set it strictly above the noise floor \eqn{\sigma\sqrt{2}} and
#'   below the smallest expected jump amplitude. When the variation of the
#'   smooth signal overlaps the jumps, select it by cross-validation.
#' @param max_iter maximum number of ADMM iterations (default 100).
#' @param rho positive augmented-Lagrangian step size; it controls the
#'   convergence speed but not the solution (default 1.0).
#'
#' @returns An object of class \code{ssj} (a list, so that all the slots below
#'   remain directly accessible) with a \code{\link[=print.ssj]{print}} and a
#'   \code{\link[=plot.ssj]{plot}} method. Its slots are:
#' \itemize{
#'   \item \code{smoothed_level}: fitted piecewise-constant signal at each
#'     (sorted) sampling site.
#'   \item \code{x}: sorted sampling sites.
#'   \item \code{nobs}: number of distinct sampling sites used.
#'   \item \code{lambda}: the \code{lambda} value used.
#'   \item \code{gamma}: absolute first differences \eqn{|u_{i+1} - u_i|} of the
#'     fit; entries greater than zero mark the detected jump locations. It has
#'     length \code{nobs - 1} and is named by the midpoints
#'     \eqn{(x_i + x_{i+1})/2}.
#' }
#'
#' @references
#' Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
#' signals. \emph{Journal of Computational and Graphical Statistics},
#' \strong{33}(2), 651--664. \doi{10.1080/10618600.2023.2262000}
#'
#' @seealso \code{\link{ssj_cssd_dp}} for the exact cubic CSSD.
#'
#' @examples
#' set.seed(1)
#' x   <- sort(runif(100))
#' y   <- sin(4 * pi * x) + c(rep(0, 50), rep(1, 50)) + rnorm(100, sd = 0.1)
#' fit <- ssj_cssd(y, x, lambda = 5)
#' plot(x, y)
#' lines(fit$x, fit$smoothed_level, col = "steelblue", lwd = 2)
#'
#' @export
ssj_cssd <- function(y, x = NULL, lambda, max_iter = 100L, rho = 1.0) {
  y <- as.numeric(y)
  n <- length(y)

  if (is.null(x)) x <- seq_len(n)
  x <- as.numeric(x)

  if (length(x) != n)       stop("'x' and 'y' must have the same length.")
  if (anyNA(y) || anyNA(x)) stop("'y' and 'x' must not contain NAs.")
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) ||
      lambda <= 0)          stop("'lambda' must be a positive scalar.")
  if (!is.numeric(rho) || length(rho) != 1L || !is.finite(rho) || rho <= 0) {
    stop("'rho' must be a positive scalar.")
  }
  max_iter <- as.integer(max_iter)
  if (is.na(max_iter) || max_iter < 1L) stop("'max_iter' must be a positive integer.")

  ord <- order(x)
  y   <- y[ord]
  x   <- x[ord]

  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    y   <- as.numeric(tapply(y, grp, mean))
    x   <- x[!duplicated(x)]
    n   <- length(x)
  }
  if (n < 2L) stop("at least two distinct values of 'x' are needed.")

  u <- as.numeric(cssd_irregular_eigen(y, x, lambda, max_iter, rho))

  jumps        <- abs(diff(u))
  names(jumps) <- (x[-n] + x[-1]) / 2

  res <- list(
    smoothed_level = u,
    x              = x,
    nobs           = n,
    lambda         = lambda,
    gamma          = jumps
  )
  new_ssj(res, method = "cssd_tv", y = y, call = match.call())
}


#' Cubic smoothing spline with discontinuities (CSSD)
#'
#' Computes the global minimiser of the weak-rod model
#' \deqn{\min_{f,J} p\sum_i\!\left(\frac{y_i-f(x_i)}{\delta_i}\right)^{\!2}
#'   + (1-p)\int(f''(t))^2\,dt + \gamma|J|}
#' using dynamic programming with incremental QR updates (Givens rotations),
#' exactly as described in Storath and Weinmann (2024).
#'
#' The algorithm computes the energy of every candidate interval, so both its
#' running time and its memory footprint grow as \eqn{O(n^2)}. Sample sizes of
#' a few thousand observations are the practical limit.
#'
#' @param y numeric vector of observations (no missing values).
#' @param x numeric vector of sampling locations (no missing values; it is
#'   sorted internally). If \code{NULL} (default), regular spacing \code{1:n} is
#'   assumed.
#' @param p stiffness parameter in \eqn{(0,1)}: weight of the data-fit term
#'   relative to the smoothness penalty. Values close to 1 favour the data fit;
#'   values close to 0 favour smoothness (approaching piecewise-linear
#'   regression as \eqn{p \to 0}). Use \code{\link{auto_ssj_cssd}} for automatic
#'   selection by K-fold cross-validation.
#' @param gamma jump-count penalty \eqn{\gamma > 0}. Larger values discourage
#'   discontinuities; \eqn{\gamma \to \infty} recovers the classical cubic
#'   smoothing spline. Use \code{\link{auto_ssj_cssd}} for automatic selection.
#' @param delta positive scalar or vector of length \code{n} with the pointwise
#'   noise standard deviations \eqn{\delta_i}. Default 1 (all equal).
#'
#' @returns An object of class \code{ssj} (a list, so that all the slots below
#'   remain directly accessible) with a \code{\link[=print.ssj]{print}} and a
#'   \code{\link[=plot.ssj]{plot}} method. Its slots are:
#' \itemize{
#'   \item \code{smoothed_level}: fitted cubic spline values at each (sorted)
#'     data site.
#'   \item \code{x}: sorted data sites.
#'   \item \code{y}: the response, sorted (and averaged over ties) to match
#'     \code{x}.
#'   \item \code{disc_locs}: x-coordinates of the detected discontinuities
#'     (midpoints between adjacent data sites).
#'   \item \code{n_disc}: number of detected discontinuities.
#'   \item \code{p}, \code{gamma}: the parameters used.
#'   \item \code{delta}: the vector of noise standard deviations used.
#'   \item \code{nobs}: number of distinct data sites used.
#' }
#'
#' @references
#' Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
#' signals. \emph{Journal of Computational and Graphical Statistics},
#' \strong{33}(2), 651--664. \doi{10.1080/10618600.2023.2262000}
#'
#' @seealso \code{\link{auto_ssj_cssd}}, \code{\link{predict_cssd_dp}},
#'   \code{\link{ssj_cssd}}.
#'
#' @examples
#' set.seed(1)
#' x   <- sort(runif(100))
#' y   <- sin(4 * pi * x) + c(rep(0, 50), rep(1, 50)) + rnorm(100, sd = 0.1)
#' fit <- ssj_cssd_dp(y, x, p = 0.999, gamma = 4)
#' plot(x, y)
#' lines(fit$x, fit$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
ssj_cssd_dp <- function(y, x = NULL, p = 0.999, gamma = 1, delta = 1) {
  y <- as.numeric(y)
  n <- length(y)

  if (is.null(x)) x <- seq_len(n)
  x <- as.numeric(x)

  if (length(x) != n)            stop("'x' and 'y' must have the same length.")
  if (anyNA(y) || anyNA(x))      stop("'y' and 'x' must not contain NAs.")
  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1) {
    stop("'p' must be a scalar in (0, 1).")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) ||
      gamma <= 0)                stop("'gamma' must be a positive scalar.")
  if (any(delta <= 0))           stop("'delta' must be positive.")

  ord <- order(x)
  y   <- y[ord]
  x   <- x[ord]

  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    y   <- as.numeric(tapply(y, grp, mean))
    x   <- x[!duplicated(x)]
    n   <- length(x)
  }
  if (n < 2L) stop("at least two distinct values of 'x' are needed.")

  delta <- rep(as.numeric(delta), length.out = n)

  raw <- cssd_dp(y, x, p = p, gamma = gamma, delta_in = delta)

  res <- list(
    smoothed_level = raw$f_hat,
    x              = x,
    y              = y,
    disc_locs      = raw$disc_locs,
    n_disc         = raw$n_disc,
    p              = p,
    gamma          = gamma,
    delta          = delta,
    nobs           = n
  )
  new_ssj(res, method = "cssd", y = y, call = match.call())
}


#' Predict from a fitted CSSD at new locations
#'
#' Evaluates the piecewise cubic spline returned by
#' \code{\link{ssj_cssd_dp}} at arbitrary query points. Query points outside
#' the range of the data are extrapolated by the cubic polynomial of the
#' nearest segment.
#'
#' @param fit object returned by \code{\link{ssj_cssd_dp}}.
#' @param xnew numeric vector of query locations (it need not be sorted).
#' @param y_orig,x_orig original data vectors, needed to re-fit the Hermite
#'   control points on each segment. If \code{NULL} (default), the sorted data
#'   stored in \code{fit} by \code{\link{ssj_cssd_dp}} are used.
#' @param delta positive scalar or vector with the noise standard deviations
#'   (the same as those used to obtain \code{fit}). If \code{NULL} (default),
#'   the value stored in \code{fit} is used.
#'
#' @returns Numeric vector of fitted values, in the same order as \code{xnew}.
#'
#' @seealso \code{\link{ssj_cssd_dp}}.
#'
#' @examples
#' set.seed(1)
#' x    <- sort(runif(100))
#' y    <- sin(4 * pi * x) + c(rep(0, 50), rep(1, 50)) + rnorm(100, sd = 0.1)
#' fit  <- ssj_cssd_dp(y, x, p = 0.999, gamma = 4)
#' xnew <- seq(0, 1, length.out = 200)
#' plot(x, y)
#' lines(xnew, predict_cssd_dp(fit, xnew), col = "darkgreen", lwd = 2)
#'
#' @export
predict_cssd_dp <- function(fit, xnew, y_orig = NULL, x_orig = NULL,
                               delta = NULL) {
  if (is.null(x_orig)) x_orig <- fit$x
  if (is.null(y_orig)) y_orig <- fit$y
  if (is.null(y_orig) || is.null(x_orig)) {
    stop("'y_orig' and 'x_orig' are needed when 'fit' does not store the data.")
  }
  n <- length(x_orig)
  if (length(y_orig) != n) stop("'x_orig' and 'y_orig' must have the same length.")
  if (is.null(delta)) delta <- if (is.null(fit$delta)) 1 else fit$delta
  delta <- rep(as.numeric(delta), length.out = n)

  ord <- order(x_orig)
  y_s <- as.numeric(y_orig)[ord]
  x_s <- as.numeric(x_orig)[ord]

  xnew   <- as.numeric(xnew)
  q_ord  <- order(xnew)
  fitted <- as.numeric(
    cssd_dp_predict(y_s, x_s, fit$smoothed_level,
                       fit$disc_locs, xnew[q_ord],
                       fit$p, delta)
  )
  # cssd_dp_predict() requires sorted query points; restore the input order
  out <- numeric(length(xnew))
  out[q_ord] <- fitted
  out
}


#' Automatic CSSD with K-fold cross-validation of the parameters
#'
#' Selects the stiffness parameter \eqn{p} and the jump-count penalty
#' \eqn{\gamma} jointly by K-fold cross-validation, as in Section 3.4 of Storath
#' and Weinmann (2024). A coarse two-dimensional grid search over
#' \eqn{(p,\gamma)} is followed by a Nelder-Mead refinement on the
#' reparameterised scale \eqn{(\mathrm{logit}(p), \mathrm{logit}(q))}, with
#' \eqn{\gamma = p\,q/(1-q)}.
#'
#' The fold assignment is drawn with \code{\link[base]{sample}}, so the result
#' depends on the state of the random number generator: call
#' \code{\link[base]{set.seed}} beforehand for reproducible output. Note also
#' that the underlying solver is \eqn{O(n^2)} in both time and memory, and that
#' the cross-validation calls it \code{k_folds} times per grid point.
#'
#' @param y numeric vector of observations (no missing values).
#' @param x numeric vector of sampling locations (no missing values). If
#'   \code{NULL} (default), regular spacing \code{1:n} is assumed.
#' @param k_folds number of cross-validation folds (default 5, as in the
#'   paper).
#' @param delta positive scalar or vector of length \code{n} with the noise
#'   standard deviations (default 1).
#' @param p_init initial stiffness (default 0.99, as in the paper).
#' @param gamma_init initial jump penalty, in the units of \code{y} (default 1).
#'   The search itself is conducted on \code{y} standardised to unit standard
#'   deviation, so that the range of jump penalties explored does not depend on
#'   the units of the response; \code{gamma_init} is standardised with it, and
#'   the \code{gamma} reported in \code{cv_par} is back in the units of
#'   \code{y}. A selected penalty at the edge of the search range draws a
#'   warning, since the number of discontinuities is then set by the range
#'   rather than by the data.
#' @param gamma_grid optional vector of candidate jump penalties for the coarse
#'   search, in units of \code{var(y)}; the default spans \code{1e-3} to
#'   \code{1e3} over 25 logarithmically spaced points, which covers the range of
#'   noise-to-signal ratios met in practice. The Nelder-Mead refinement is not
#'   confined to it, but a coarse grid that misses the optimum by orders of
#'   magnitude can leave the refinement in a local dip.
#' @param maxit_nm maximum number of Nelder-Mead iterations (default 200).
#'
#' @returns The output of \code{\link{ssj_cssd_dp}} at the cross-validated
#'   optimal \eqn{(p,\gamma)}, with two extra slots:
#' \itemize{
#'   \item \code{cv_score}: the minimised cross-validation score.
#'   \item \code{cv_par}: named vector with the selected \code{p} and
#'     \code{gamma}.
#' }
#'
#' @references
#' Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
#' signals. \emph{Journal of Computational and Graphical Statistics},
#' \strong{33}(2), 651--664. \doi{10.1080/10618600.2023.2262000}
#'
#' @seealso \code{\link{ssj_cssd_dp}}, \code{\link{predict_cssd_dp}}.
#'
#' @examples
#' set.seed(1)
#' x   <- sort(runif(60))
#' y   <- sin(4 * pi * x) + c(rep(0, 30), rep(1, 30)) + rnorm(60, sd = 0.1)
#' fit <- auto_ssj_cssd(y, x)
#' fit$cv_par
#' plot(x, y)
#' lines(fit$x, fit$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
auto_ssj_cssd <- function(y, x = NULL, k_folds = 5L, delta = 1,
                          p_init = 0.99, gamma_init = 1, gamma_grid = NULL,
                          maxit_nm = 200L) {
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(x)) x <- seq_len(n)
  x <- as.numeric(x)
  if (length(x) != n)       stop("'x' and 'y' must have the same length.")
  if (anyNA(y) || anyNA(x)) stop("'y' and 'x' must not contain NAs.")

  ord <- order(x)
  y   <- y[ord]
  x   <- x[ord]
  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    y   <- as.numeric(tapply(y, grp, mean))
    x   <- x[!duplicated(x)]
    n   <- length(x)
  }
  if (n < 2L) stop("at least two distinct values of 'x' are needed.")
  delta <- rep(as.numeric(delta), length.out = n)

  # The search is conducted on y standardised to unit standard deviation, and
  # the selected gamma is mapped back at the end.  gamma is a charge per
  # discontinuity levied against the data-fit and smoothness terms, both of
  # which scale with the square of the units of y, so a grid fixed in absolute
  # terms is meaningful for exactly one scaling of the response.  Searching
  # [0.5, 100] with y in MWh (sd ~ 1.5e3, so sigma^2 ~ 1.3e6) puts the whole
  # grid at under 1e-4 of a single squared residual, where every candidate buys
  # a discontinuity for nothing: on the electricity-demand data that produced
  # 146 discontinuities among 295 sites, at a cross-validated error 2.8 times
  # worse than the optimum the search never reached.  The same grid on a signal
  # of unit scale is the opposite extreme.  Both terms of the objective carry
  # the factor, so p is unaffected and only gamma has to be rescaled.
  y_scale <- stats::sd(y)
  if (!is.finite(y_scale) || y_scale <= 0) y_scale <- 1
  yz <- y / y_scale

  fold_id <- sample(rep(seq_len(k_folds), length.out = n))

  # Paper's reparameterisation (Section 3.4): optimise over (p, q) where
  # gamma = p * q / (1 - q).  This mirrors dividing the objective by p and
  # mapping gamma' = gamma/p via q/(1-q), giving a compact search space.
  # par = [logit(p), logit(q)]
  # Clamp q away from 1 as well as p: gamma = p*q/(1-q) needs logit(q) ~ 35 to
  # reach gamma ~ 1e13, and beyond that 1 - q underflows to zero and gamma
  # becomes Inf, which the dynamic program cannot use.  Standardising the
  # response keeps logit(q) near 2.5 whatever the units of y, so the optimiser
  # never goes near this on its own; the guard is for a 'gamma_grid' supplied
  # in the extreme.  q = 1 - 1e-12 already means a penalty of order 1e12 times
  # var(y), which is the classical cubic spline for any practical purpose.
  q_max <- 1 - 1e-12
  cv_score <- function(par) {
    p_try <- 1 / (1 + exp(-par[1]))
    q_try <- min(1 / (1 + exp(-par[2])), q_max)
    # Clamp p away from 1: beta = sqrt(1-p) -> 0 causes QR degeneracy
    p_try   <- min(p_try, 0.9999)
    gam_try <- p_try * q_try / (1 - q_try)
    err <- 0
    for (fold in seq_len(k_folds)) {
      tr <- which(fold_id != fold)
      te <- which(fold_id == fold)
      if (length(tr) < 3L || length(te) == 0L) next
      fit_tr <- tryCatch(
        cssd_dp(yz[tr], x[tr], p = p_try, gamma = gam_try,
                   delta_in = delta[tr]),
        error = function(e) NULL
      )
      if (is.null(fit_tr)) return(1e12)
      pred <- as.numeric(
        cssd_dp_predict(yz[tr], x[tr], fit_tr$f_hat,
                           fit_tr$disc_locs, x[te],
                           p_try, delta[tr])
      )
      err <- err + sum(((yz[te] - pred) / delta[te])^2, na.rm = TRUE)
    }
    err / n
  }

  # Coarse 2-D grid search (p x gamma), then Nelder-Mead refinement.
  # This replaces Matlab's simulannealbnd: bounded, deterministic, and equally
  # reliable for this problem size (Section 3.4 notes grid search as the
  # parallel-processing alternative to simulated annealing).
  # gam_grid is in the units of the STANDARDISED response, and spans six orders
  # of magnitude: the optimum moves with the noise-to-signal ratio, and a narrow
  # grid whose end is binding is silently satisfied rather than reported.  On the
  # Section-3 design at sigma = 0.05 the cross-validated optimum sits at 0.33 in
  # these units, below the old floor of 0.5, which is why that estimator looked
  # more conservative there than cross-validation wanted.
  p_grid   <- c(0.9, 0.99, 0.999)
  gam_grid <- if (is.null(gamma_grid))
    exp(seq(log(1e-3), log(1e3), length.out = 25L))
  else sort(as.numeric(gamma_grid))
  if (!length(gam_grid) || any(!is.finite(gam_grid)) || any(gam_grid <= 0))
    stop("'gamma_grid' must be a vector of positive, finite values")

  # Initial (p, q) from (p_init, gamma_init): q = gamma / (p + gamma).
  # gamma_init is given in the units of y, so it is standardised like the grid.
  gamma_init <- gamma_init / y_scale^2
  q_init   <- gamma_init / (p_init + gamma_init)
  best_val <- Inf
  best_par <- c(log(p_init / (1 - p_init)), log(q_init / (1 - q_init)))
  for (pv in p_grid) {
    for (gv in gam_grid) {
      qv  <- gv / (pv + gv)
      par <- c(log(pv / (1 - pv)), log(qv / (1 - qv)))
      val <- cv_score(par)
      if (val < best_val) { best_val <- val; best_par <- par }
    }
  }

  # Local refinement from the best grid point
  opt <- stats::optim(best_par, cv_score, method = "Nelder-Mead",
                      control = list(maxit = maxit_nm, reltol = 1e-4))

  p_opt   <- min(1 / (1 + exp(-opt$par[1])), 0.9999)
  q_opt   <- min(1 / (1 + exp(-opt$par[2])), q_max)
  gam_z   <- p_opt * q_opt / (1 - q_opt)      # standardised units
  gam_opt <- gam_z * y_scale^2                # back to the units of y

  # A selected gamma at either end of the grid means the grid, not the data, is
  # deciding how many discontinuities there are.  Say so rather than returning a
  # constrained answer silently -- that is how the electricity-demand fit came
  # to report 146 discontinuities without complaint.
  if (gam_z <= 1.5 * min(gam_grid) || gam_z >= max(gam_grid) / 1.5)
    warning("auto_ssj_cssd() selected a jump penalty at the edge of its search ",
            "range (gamma = ", signif(gam_opt, 4), ", i.e. ", signif(gam_z, 4),
            " in units of var(y)): the number of discontinuities may be set by ",
            "the search range rather than by the data.", call. = FALSE)

  fit <- ssj_cssd_dp(y, x, p = p_opt, gamma = gam_opt, delta = delta)
  # cv_score is reported in the units of y^2, like the data it came from
  fit$cv_score <- opt$value * y_scale^2
  fit$cv_par   <- c(p = p_opt, gamma = gam_opt)
  fit$call     <- match.call()
  fit
}
