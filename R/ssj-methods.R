# S3 class for the smoothing splines with jumps returned by the ssj_* family,
# together with its print and plot methods.


# Human-readable labels for the estimation methods.  Kept in one place so that
# print() and plot() cannot drift apart.  Note that "cssd" and "cssd_tv" are
# NOT the estimators of Pelagatti, Maranzano and Gatto Ronchieri: they are the
# competing proposal of Storath and Weinmann (2024), included for comparison,
# and the labels say so explicitly.
ssj_method_label <- function(method) {
  switch(method,
    mle     = "SSJ-MLE: domain warping, state-space maximum likelihood",
    kfs     = "SSJ-KFS: domain warping, Kalman filter and smoother, lambda fixed",
    atw     = "SSJ-ATW: alternating domain-warping algorithm",
    fix     = "SSJ-FIX: time-varying disturbance variances, lambda fixed",
    cssd    = "CSSD: cubic smoothing spline with discontinuities (Storath and Weinmann, 2024)",
    cssd_tv = "TV denoising: piecewise-constant approximation to CSSD",
    method
  )
}

# Is this one of the domain-warping estimators of the paper, i.e. does the
# object carry a warping vector gamma of length n-1?
ssj_is_warp <- function(object) {
  object$method %in% c("mle", "kfs", "atw")
}

# Internal constructor.  It never drops or renames anything already present in
# `out`, so code written against the plain lists returned by earlier versions
# keeps working; it only adds the class and a few derived fields.
new_ssj <- function(out, method, y, call) {
  out$method <- method
  out$y      <- y
  out$call   <- call
  structure(out, class = "ssj")
}


#' Discontinuities declared by a smoothing spline with jumps
#'
#' Applies the detection rule of Pelagatti, Maranzano and Gatto Ronchieri
#' (2026) to a fitted object of class \code{ssj} and returns the declared
#' discontinuities, their locations, the corresponding domain warps and the
#' size of the jump in the fitted curve.
#'
#' For the domain-warping estimators (\code{\link{ssj_mle}},
#' \code{\link{ssj_kfs}}, \code{\link{ssj_atw}}) an interval is declared a
#' discontinuity when the estimated warp \eqn{\hat\gamma_t} clears \emph{both}
#' \code{frac_total} times the total warping in the sample and
#' \code{frac_delta} times the mean sampling interval. The first condition
#' filters out numerically small, spurious warps, the second keeps the rule
#' meaningful in the units of \code{x}; together they retain genuine jumps of
#' very different magnitudes. Supplying \code{thresh} replaces both conditions
#' by the single rule \eqn{\hat\gamma_t > \code{thresh}}, which is convenient
#' with equally spaced data, where one unit of the domain has a direct
#' interpretation (a week, say) and \code{thresh = 1} means that the algorithm
#' has inserted at least one extra unit of fictitious domain length.
#'
#' Since \eqn{\gamma_t} lives on the interval between \eqn{x_t} and
#' \eqn{x_{t+1}}, a break is dated at \eqn{x_{t+1}}, the first observation
#' belonging to the new regime.
#'
#' For \code{\link{ssj_fix}}, which parameterises jumps through time-varying
#' disturbance standard deviations rather than domain warps, only the
#' \code{frac_total} rule is applied, because the \code{frac_delta} condition
#' has no meaning in those units. For the CSSD estimators the discontinuity set
#' is part of the solution and is returned as it is, ignoring the thresholds.
#'
#' @param object an object of class \code{ssj}.
#' @param frac_total minimum share of the total warping that an interval must
#'   carry to be declared a discontinuity (default 0.05).
#' @param frac_delta minimum multiple of the mean sampling interval that the
#'   warp must reach (default 0.30).
#' @param thresh optional absolute threshold on the warp; when supplied it
#'   overrides \code{frac_total} and \code{frac_delta}.
#'
#' @returns A \code{data.frame}, with one row per declared discontinuity, and
#'   columns:
#' \itemize{
#'   \item \code{index}: index of the first observation of the new regime.
#'   \item \code{x}: its location, that is the date of the break.
#'   \item \code{warp}: the estimated domain warp \eqn{\hat\gamma_t} (or the
#'     additional standard deviation for \code{\link{ssj_fix}}, \code{NA} for
#'     the CSSD estimators).
#'   \item \code{jump}: the size of the jump, that is the change of the fitted
#'     curve across the break.
#' }
#' It has zero rows when no discontinuity is declared.
#'
#' @references
#' Pelagatti, M., Maranzano, P. and Gatto Ronchieri, G. (2026). Smoothing
#' splines with some discontinuities.
#'
#' @seealso \code{\link{plot.ssj}}, \code{\link{print.ssj}}.
#'
#' @examples
#' # 'maxsum' is a warping budget in the units of x; the range of eruptions
#' # is about 3.5, so a budget of 1 is already generous.
#' fit <- ssj_atw(y = faithful$waiting, x = faithful$eruptions,
#'                lambda = 1, maxsum = 1)
#' ssj_breaks(fit)
#'
#' @export
ssj_breaks <- function(object, frac_total = 0.05, frac_delta = 0.30,
                       thresh = NULL) {
  if (!inherits(object, "ssj")) stop("'object' must be of class \"ssj\"")

  xx  <- object$x
  fit <- object$smoothed_level
  n   <- length(xx)

  empty <- data.frame(index = integer(0), x = numeric(0),
                      warp = numeric(0), jump = numeric(0))

  # CSSD (Storath and Weinmann): the discontinuity set is part of the solution
  if (identical(object$method, "cssd")) {
    locs <- object$disc_locs
    if (is.null(locs) || length(locs) == 0L) return(empty)
    idx <- vapply(locs, function(l) which(xx > l)[1L], integer(1))
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) return(empty)
    return(data.frame(index = idx, x = xx[idx], warp = NA_real_,
                      jump = fit[idx] - fit[idx - 1L]))
  }

  # TV denoising: the ADMM iterate is only approximately sparse, so a bare
  # numerical tolerance would flag almost every interval.  Apply the same
  # share-of-total rule used for the warps (frac_delta does not apply here,
  # because these differences are in the units of y, not of x).
  if (identical(object$method, "cssd_tv")) {
    g <- as.numeric(object$gamma)
    g[!is.finite(g)] <- 0
    sel <- if (!is.null(thresh)) which(g > thresh) else which(g > frac_total * sum(g))
    if (length(sel) == 0L) return(empty)
    idx <- sel + 1L
    return(data.frame(index = idx, x = xx[idx], warp = g[sel],
                      jump = fit[idx] - fit[idx - 1L]))
  }

  if (ssj_is_warp(object)) {
    g <- as.numeric(object$gamma)          # length n-1, warp on (x_t, x_{t+1})
  } else {                                  # ssj_fix: additional std devs
    g <- as.numeric(object$sigmas)[-1L]     # first element is NA by construction
  }
  if (is.null(g) || length(g) == 0L || all(!is.finite(g))) return(empty)
  g[!is.finite(g)] <- 0

  if (!is.null(thresh)) {
    sel <- which(g > thresh)
  } else {
    tot <- sum(g)
    ok  <- g > frac_total * tot
    if (ssj_is_warp(object)) {
      # the second condition is expressed in the units of x, so it only makes
      # sense for the domain-warping estimators
      ok <- ok & (g > frac_delta * mean(diff(xx)))
    }
    sel <- which(ok)
  }
  if (length(sel) == 0L) return(empty)

  idx <- sel + 1L                           # first observation of the new regime
  idx[idx > n] <- n
  data.frame(index = idx, x = xx[idx], warp = g[sel],
             jump = fit[idx] - fit[pmax(idx - 1L, 1L)])
}


#' Print method for the class ssj
#'
#' Prints the estimation method, the call, the estimated parameters, the
#' information criteria and a table of the declared discontinuities.
#'
#' @param x an object of class \code{ssj}.
#' @param digits number of significant digits used for the numbers printed.
#' @param ... further arguments passed to \code{\link{ssj_breaks}}, such as
#'   \code{thresh}, \code{frac_total} or \code{frac_delta}.
#'
#' @returns \code{x}, invisibly. Called for its side effect of printing.
#'
#' @seealso \code{\link{plot.ssj}}, \code{\link{ssj_breaks}}.
#'
#' @examples
#' # 'maxsum' is a warping budget in the units of x; the range of eruptions
#' # is about 3.5, so a budget of 1 is already generous.
#' fit <- ssj_atw(y = faithful$waiting, x = faithful$eruptions,
#'                lambda = 1, maxsum = 1)
#' print(fit)
#'
#' @export
print.ssj <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Smoothing spline with jumps\n")
  cat("Method: ", ssj_method_label(x$method), "\n", sep = "")
  if (!is.null(x$call)) {
    cat("Call:\n  ", paste(deparse(x$call), collapse = "\n   "), "\n", sep = "")
  }

  # NB: exact name tests throughout.  `x$p` would partial-match `pars`, because
  # `$` on a list falls back to partial matching when there is no exact hit.
  has <- function(nm) nm %in% names(x)
  say <- function(nm, value, note = "")
    cat("  ", format(nm, width = 8), " = ", format(value, digits = digits),
        note, "\n", sep = "")

  cat("\nParameters:\n")
  if (has("pars")) for (nm in names(x$pars)) say(nm, x$pars[[nm]])
  if (has("p"))                                        say("p", x[["p"]])
  if (identical(x$method, "cssd")    && has("gamma"))  say("gamma",  x[["gamma"]])
  if (identical(x$method, "cssd_tv") && has("lambda")) say("lambda", x[["lambda"]])
  if (has("maxsum")) {
    note <- if (identical(x$method, "fix"))
      "   (maximum sum of the additional standard deviations)"
    else
      "   (maximum sum of the domain warps)"
    say("budget", x[["maxsum"]], note)
  }

  cat("\nObservations: ", x$nobs, "\n", sep = "")
  if (!is.null(x$df))
    cat("Effective degrees of freedom: ", format(x$df, digits = digits), "\n", sep = "")
  if (!is.null(x$loglik))
    cat("Log-likelihood: ", format(x$loglik, digits = digits), "\n", sep = "")
  if (!is.null(x$ic)) {
    cat("Information criteria:\n")
    nms <- names(x$ic)
    for (i in seq_along(x$ic))
      cat("  ", format(toupper(nms[i]), width = 5), " = ",
          format(x$ic[[i]], digits = digits), "\n", sep = "")
  }

  br <- ssj_breaks(x, ...)
  cat("\nDeclared discontinuities: ", nrow(br), "\n", sep = "")
  if (nrow(br) > 0L) {
    tab <- data.frame(x = format(br$x, digits = digits),
                      jump = format(br$jump, digits = digits))
    if (any(is.finite(br$warp))) tab$warp <- format(br$warp, digits = digits)
    row.names(tab) <- format(br$index)
    print(tab, quote = FALSE)
  }
  invisible(x)
}


#' Plot method for the class ssj
#'
#' Draws the observed data with the fitted smoothing spline with jumps and,
#' optionally, a confidence band and vertical lines at the declared
#' discontinuities. For the domain-warping estimators a second panel below
#' shows the estimated warps \eqn{\hat\gamma_t}, which measure the strength of
#' each jump, together with the detection threshold. This reproduces Figure 4
#' of Pelagatti, Maranzano and Gatto Ronchieri (2026).
#'
#' @param x an object of class \code{ssj}.
#' @param y not used: for consistency with the generic function.
#' @param prob coverage probability of the confidence band around the fitted
#'   curve; if \code{NULL} (default), or if the estimator does not return
#'   pointwise variances, no band is drawn.
#' @param show_breaks logical: if \code{TRUE} (default), vertical dashed lines
#'   are drawn at the declared discontinuities.
#' @param show_warp logical: if \code{TRUE} (default), a lower panel with the
#'   estimated domain warps is added, when the estimator provides them.
#' @param main title of the plot; the default, \code{NULL}, draws no title.
#' @param xlab,ylab axis labels; if \code{NULL} defaults are used.
#' @param col_fit,col_data colours of the fitted curve and of the observations.
#' @param ... further arguments passed to \code{\link{ssj_breaks}}, such as
#'   \code{thresh}, \code{frac_total} or \code{frac_delta}.
#'
#' @returns \code{NULL}, invisibly. Called for its side effect of drawing.
#'
#' @references
#' Pelagatti, M., Maranzano, P. and Gatto Ronchieri, G. (2026). Smoothing
#' splines with some discontinuities.
#'
#' @seealso \code{\link{print.ssj}}, \code{\link{ssj_breaks}}.
#'
#' @examples
#' # 'maxsum' is a warping budget in the units of x; the range of eruptions
#' # is about 3.5, so a budget of 1 is already generous.
#' fit <- ssj_atw(y = faithful$waiting, x = faithful$eruptions,
#'                lambda = 1, maxsum = 1)
#' plot(fit)
#'
#' # no title is drawn unless one is asked for
#' plot(fit, main = "Old Faithful", xlab = "eruption length",
#'      ylab = "waiting time")
#'
#' # a single absolute threshold on the warps, as in Section 4.2 of the paper
#' plot(fit, thresh = 1)
#'
#' @export
plot.ssj <- function(x, y = NULL, prob = NULL, show_breaks = TRUE,
                     show_warp = TRUE, main = NULL, xlab = NULL, ylab = NULL,
                     col_fit = "firebrick", col_data = "steelblue", ...) {
  obj <- x
  xx  <- obj$x
  yy  <- obj$y
  fit <- obj$smoothed_level

  br      <- ssj_breaks(obj, ...)
  has_warp <- (ssj_is_warp(obj) || identical(obj$method, "fix")) &&
              !is.null(if (ssj_is_warp(obj)) obj$gamma else obj$sigmas)
  panel2  <- isTRUE(show_warp) && has_warp

  # No title by default: the method labels are long, and a plot that is going
  # into a paper usually wants either no title or one the user chooses.
  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "y"

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  if (panel2) {
    # reserve room at the top only when there is actually a title to draw
    par(mfrow = c(2, 1), mar = c(0.5, 4.5, if (is.null(main)) 1 else 3, 1),
        oma = c(4, 0, 0, 0))
  }

  # ---- upper panel: data, fitted curve, declared discontinuities ----------
  band <- NULL
  if (!is.null(prob) && !is.null(obj$var_smoothed_level)) {
    z    <- abs(qnorm((1 - prob)/2))
    se   <- sqrt(pmax(obj$var_smoothed_level, 0))
    band <- list(lo = fit - z*se, up = fit + z*se)
  }
  ylim <- range(c(yy, fit, band$lo, band$up), na.rm = TRUE)

  plot(xx, if (is.null(yy)) fit else yy, type = "n", ylim = ylim,
       xaxt = if (panel2) "n" else "s",
       xlab = if (panel2) "" else xlab, ylab = ylab, main = main)
  if (!is.null(band)) {
    polygon(c(xx, rev(xx)), c(band$up, rev(band$lo)),
            col = adjustcolor(col_fit, alpha.f = 0.20), border = NA)
  }
  if (!is.null(yy)) {
    points(xx, yy, pch = 16, cex = 0.5,
           col = adjustcolor(col_data, alpha.f = 0.35))
  }
  if (isTRUE(show_breaks) && nrow(br) > 0L)
    abline(v = br$x, lty = 2, col = "grey40")
  lines(xx, fit, col = col_fit, lwd = 2)

  # ---- lower panel: strength of the jumps --------------------------------
  if (panel2) {
    g   <- if (ssj_is_warp(obj)) as.numeric(obj$gamma) else as.numeric(obj$sigmas)[-1L]
    g[!is.finite(g)] <- 0
    xg  <- xx[-1L]
    lab <- if (ssj_is_warp(obj)) expression("Domain warp" ~ gamma) else "Additional s.d."

    par(mar = c(0.5, 4.5, 0.5, 1))
    plot(xg, g, type = "h", col = "grey80", xlab = "", ylab = lab)
    if (isTRUE(show_breaks) && nrow(br) > 0L) {
      abline(v = br$x, lty = 2, col = "grey40")
      lines(br$x, br$warp, type = "h", col = col_fit, lwd = 2)
      points(br$x, br$warp, pch = 16, col = col_fit, cex = 0.7)
      cut <- min(br$warp)
      if (is.finite(cut)) abline(h = cut, lty = 3, col = "grey40")
    }
    mtext(xlab, side = 1, outer = TRUE, line = 2.5)
  }

  invisible(NULL)
}
