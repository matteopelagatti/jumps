#' TV denoising (piecewise-constant smoothing) for irregularly sampled signals
#'
#' Solves the total-variation denoising problem
#' \deqn{\min_u \frac{1}{2}\|y - u\|_2^2 + \lambda\|Du\|_1,}
#' where \eqn{D} is the unweighted first-difference operator
#' \eqn{(Du)_i = u_{i+1} - u_i}, using ADMM with a pre-factorized
#' \eqn{I + \rho D^\top D} linear system.  The minimizer \eqn{u} is piecewise
#' constant, with nonzero first differences at detected jump locations.
#'
#' This implements the simplified TV / first-order variant related to the CSSD
#' framework of Storath and Weinmann (2024). The full CSSD (piecewise cubic
#' smoothing spline with discontinuities) requires a specialized QR-update
#' solver; the TV formulation here replaces the \eqn{C^2} roughness penalty
#' with an \eqn{\ell_1} penalty on first differences.
#'
#' The unweighted \eqn{D} (entries \eqn{\pm 1}) is used rather than the
#' slope-based \eqn{(Du)_i = (u_{i+1}-u_i)/\delta_i}: the \eqn{1/\delta}
#' weighting inflates the eigenvalues of \eqn{D^\top D} by \eqn{1/h^2 \sim
#' 10^4} for \eqn{h\sim 0.01}, which causes the ADMM \eqn{u}-update to drive
#' \eqn{u} to nearly zero in the first iteration and stall permanently.
#'
#' @param y numeric vector of observations (no missing values).
#' @param x numeric vector of sampling locations (no missing values; need
#'   not be sorted). If \code{NULL} (default), assumes regular spacing
#'   \code{1:n}.
#' @param lambda positive regularization parameter controlling the
#'   soft-threshold \eqn{\lambda/\rho} applied to first differences of
#'   \eqn{u}.  Set it strictly above the noise floor
#'   \eqn{\sigma\sqrt{2}} and below the smallest expected jump amplitude.
#'   Use \code{\link{cssd_cv_lambda}} or 5-fold cross-validation for
#'   automatic selection when the smooth-signal variation overlaps the jumps.
#' @param max_iter maximum number of ADMM iterations (default 100).
#' @param rho positive augmented-Lagrangian step size; controls convergence
#'   speed but does not affect the solution (default 1.0).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{smoothed_level}: fitted piecewise-constant signal at each
#'     (sorted) sample site.
#'   \item \code{x}: sorted sample sites.
#'   \item \code{nobs}: number of observations used.
#'   \item \code{lambda}: the \code{lambda} value used.
#'   \item \code{gamma}: absolute first differences \eqn{|u_{i+1} - u_i|}
#'     of the fit; entries greater than zero mark detected jump locations.
#'     Length \code{nobs - 1}, named by the midpoints \eqn{(x_i + x_{i+1})/2}.
#' }
#'
#' @references
#' Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
#' signals. \emph{Journal of Computational and Graphical Statistics},
#' \strong{33}(2), 651--664. \doi{10.1080/10618600.2023.2262000}
#'
#' @examples
#' x   <- sort(runif(100))
#' y   <- sin(4 * pi * x) + c(rep(0, 50), rep(1, 50)) + rnorm(100, sd = 0.1)
#' fit <- ssj_cssd(y, x, lambda = 5)
#' plot(x, y); lines(fit$x, fit$smoothed_level, col = "steelblue", lwd = 2)
#'
#' @export
ssj_cssd <- function(y, x = NULL, lambda, max_iter = 100L, rho = 1.0) {
  y <- as.numeric(y)
  n <- length(y)

  if (is.null(x)) x <- seq_len(n)
  x <- as.numeric(x)

  if (anyNA(y) || anyNA(x)) stop("y and x must not contain NAs.")
  if (length(x) != n)       stop("x and y must have the same length.")

  ord <- order(x)
  y   <- y[ord]
  x   <- x[ord]

  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    y   <- as.numeric(tapply(y, grp, mean))
    x   <- x[!duplicated(x)]
    n   <- length(x)
  }

  u <- as.numeric(cssd_irregular_eigen(y, x, lambda, as.integer(max_iter), rho))

  jumps        <- abs(diff(u))
  names(jumps) <- (x[-n] + x[-1]) / 2

  list(
    smoothed_level = u,
    x              = x,
    nobs           = n,
    lambda         = lambda,
    gamma          = jumps
  )
}


#' Cubic smoothing spline with discontinuities (paper's exact CSSD)
#'
#' Computes the global minimiser of the weak-rod model
#' \deqn{\min_{f,J} p\sum_i\!\left(\frac{y_i-f(x_i)}{\delta_i}\right)^{\!2}
#'   + (1-p)\int(f''(t))^2\,dt + \gamma|J|}
#' using dynamic programming with incremental QR updates (Givens rotations),
#' exactly as described in Storath & Weinmann (2024).
#'
#' @param y numeric vector of observations (no NAs).
#' @param x numeric vector of sampling locations (no NAs; sorted internally).
#' @param p stiffness parameter in \eqn{(0,1)}: weight of the data-fit term
#'   relative to the smoothness penalty. Values close to 1 favour data fit;
#'   values close to 0 favour smoothness (approaching piecewise-linear
#'   regression as \eqn{p \to 0}). Use \code{\link{auto_ssj_cssd}} for
#'   automatic selection via K-fold CV.
#' @param gamma jump-count penalty \eqn{\gamma > 0}. Larger values discourage
#'   discontinuities; \eqn{\gamma = \infty} recovers the classical cubic
#'   smoothing spline. Use \code{\link{auto_ssj_cssd}} for automatic selection.
#' @param delta positive scalar or vector of length \code{n}: pointwise noise
#'   standard deviations \eqn{\delta_i}. Default 1 (all equal).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{smoothed_level}: fitted cubic spline values at each (sorted)
#'     data site.
#'   \item \code{x}: sorted data sites.
#'   \item \code{disc_locs}: x-coordinates of detected discontinuities
#'     (midpoints between adjacent data sites).
#'   \item \code{n_disc}: number of detected discontinuities.
#'   \item \code{p}, \code{gamma}: parameters used.
#'   \item \code{nobs}: number of observations.
#' }
#'
#' @references
#' Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
#' signals. \emph{Journal of Computational and Graphical Statistics},
#' \strong{33}(2), 651--664. \doi{10.1080/10618600.2023.2262000}
#'
#' @examples
#' x   <- sort(runif(100))
#' y   <- sin(4 * pi * x) + c(rep(0, 50), rep(1, 50)) + rnorm(100, sd = 0.1)
#' fit <- ssj_cssd_potts(y, x, p = 0.999, gamma = 4)
#' plot(x, y); lines(fit$x, fit$smoothed_level, col = "darkgreen", lwd = 2)
#'
#' @export
ssj_cssd_potts <- function(y, x = NULL, p = 0.999, gamma = 1, delta = 1) {
  y <- as.numeric(y)
  n <- length(y)

  if (is.null(x)) x <- seq_len(n)
  x <- as.numeric(x)

  if (anyNA(y) || anyNA(x))      stop("y and x must not contain NAs.")
  if (length(x) != n)            stop("x and y must have the same length.")
  if (p <= 0 || p >= 1)          stop("p must be in (0, 1).")
  if (gamma <= 0)                 stop("gamma must be positive.")

  ord <- order(x)
  y   <- y[ord]
  x   <- x[ord]

  if (anyDuplicated(x)) {
    grp <- cumsum(c(TRUE, diff(x) != 0))
    y   <- as.numeric(tapply(y, grp, mean))
    x   <- x[!duplicated(x)]
    n   <- length(x)
  }

  delta <- rep(as.numeric(delta), length.out = n)

  raw <- cssd_potts(y, x, p = p, gamma = gamma, delta_in = delta)

  list(
    smoothed_level = raw$f_hat,
    x              = x,
    disc_locs      = raw$disc_locs,
    n_disc         = raw$n_disc,
    p              = p,
    gamma          = gamma,
    nobs           = n
  )
}


#' Predict from a fitted CSSD on new x locations
#'
#' Evaluates the piecewise cubic spline returned by \code{\link{ssj_cssd_potts}}
#' at arbitrary query points.
#'
#' @param fit   object returned by \code{ssj_cssd_potts}.
#' @param xnew  numeric vector of query locations.
#' @param y_orig,x_orig  original data vectors (needed for re-fitting Hermite
#'   control points on each segment).
#' @param delta positive scalar or vector: noise std devs (same as used to
#'   fit \code{fit}).
#'
#' @return Numeric vector of fitted values at \code{xnew}.
#' @export
predict_cssd_potts <- function(fit, xnew, y_orig, x_orig, delta = 1) {
  n     <- length(x_orig)
  delta <- rep(as.numeric(delta), length.out = n)
  ord   <- order(x_orig)
  y_s   <- y_orig[ord]
  x_s   <- x_orig[ord]
  as.numeric(
    cssd_potts_predict(y_s, x_s, fit$smoothed_level,
                        fit$disc_locs, sort(xnew),
                        fit$p, delta)
  )
}


#' Automatic CSSD with K-fold cross-validation parameter selection
#'
#' Selects the stiffness parameter \eqn{p} and jump-count penalty \eqn{\gamma}
#' jointly by K-fold cross-validation, as in Section 3.4 of Storath &
#' Weinmann (2024). Optimisation uses Nelder-Mead over
#' \eqn{(\log p/(1-p),\, \log\gamma)} starting from \eqn{p_0 = 0.99,
#' \gamma_0 = 1}.
#'
#' @param y numeric vector of observations (no NAs).
#' @param x numeric vector of sampling locations.
#' @param k_folds number of CV folds (default 5, as in the paper).
#' @param delta noise std dev: scalar or vector of length \code{n} (default 1).
#' @param p_init initial stiffness (default 0.99, as in the paper).
#' @param gamma_init initial jump penalty (default 1, as in the paper).
#' @param maxit maximum Nelder-Mead iterations (default 200).
#'
#' @return The output of \code{ssj_cssd_potts} at the CV-optimal \eqn{(p,\gamma)}.
#' @export
auto_ssj_cssd <- function(y, x, k_folds = 5L, delta = 1,
                           p_init = 0.99, gamma_init = 1,
                           maxit_nm = 200L) {
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(x)) x <- seq_len(n)
  x     <- as.numeric(x)
  ord   <- order(x)
  y     <- y[ord]
  x     <- x[ord]
  delta <- rep(as.numeric(delta), length.out = n)

  fold_id <- sample(rep(seq_len(k_folds), length.out = n))

  # Paper's reparameterisation (Section 3.4): optimise over (p, q) where
  # gamma = p * q / (1 - q).  This mirrors dividing the objective by p and
  # mapping gamma' = gamma/p via q/(1-q), giving a compact search space.
  # par = [logit(p), logit(q)]
  cv_score <- function(par) {
    p_try <- 1 / (1 + exp(-par[1]))
    q_try <- 1 / (1 + exp(-par[2]))
    # Clamp p away from 1: beta = sqrt(1-p) -> 0 causes QR degeneracy
    p_try   <- min(p_try, 0.9999)
    gam_try <- p_try * q_try / (1 - q_try)
    err <- 0
    for (fold in seq_len(k_folds)) {
      tr <- which(fold_id != fold)
      te <- which(fold_id == fold)
      if (length(tr) < 3L) next
      fit_tr <- tryCatch(
        cssd_potts(y[tr], x[tr], p = p_try, gamma = gam_try,
                   delta_in = delta[tr]),
        error = function(e) NULL
      )
      if (is.null(fit_tr)) return(1e12)
      pred <- as.numeric(
        cssd_potts_predict(y[tr], x[tr], fit_tr$f_hat,
                           fit_tr$disc_locs, x[te],
                           p_try, delta[tr])
      )
      err <- err + sum(((y[te] - pred) / delta[te])^2, na.rm = TRUE)
    }
    err / n
  }

  # Initial (p, q) from (p_init, gamma_init): q = gamma / (p + gamma)
  # Coarse 2-D grid search (p × gamma), then Nelder-Mead refinement.
  # This replaces Matlab's simulannealbnd: bounded, deterministic, equally
  # reliable for this problem size (Section 3.4 notes grid search as the
  # parallel-processing alternative to simulated annealing).
  p_grid   <- c(0.9, 0.99, 0.999)
  gam_grid <- exp(seq(log(0.5), log(100), length.out = 20L))

  best_val <- Inf
  best_par <- c(log(p_init / (1 - p_init)),
                log(gamma_init / (p_init * (1 + gamma_init / p_init) - gamma_init)))
  for (pv in p_grid) {
    for (gv in gam_grid) {
      qv  <- gv / (pv + gv)
      par <- c(log(pv / (1 - pv)), log(qv / (1 - qv)))
      val <- cv_score(par)
      if (val < best_val) { best_val <- val; best_par <- par }
    }
  }

  # Local refinement from best grid point
  opt <- stats::optim(best_par, cv_score, method = "Nelder-Mead",
                      control = list(maxit = maxit_nm, reltol = 1e-4))

  p_opt   <- min(1 / (1 + exp(-opt$par[1])), 0.9999)
  q_opt   <- 1 / (1 + exp(-opt$par[2]))
  gam_opt <- p_opt * q_opt / (1 - q_opt)

  fit <- ssj_cssd_potts(y, x, p = p_opt, gamma = gam_opt, delta = delta)
  fit$cv_score <- opt$value
  fit$cv_par   <- c(p = p_opt, gamma = gam_opt)
  fit
}
