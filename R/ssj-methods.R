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

# Validate and recycle observation weights.  Returns a numeric vector of
# length n; NULL means unit weights.  Weights are interpreted as counts of
# independent observations, so the observation variance is divided by them.
ssj_check_weights <- function(weights, n) {
  if (is.null(weights)) return(rep(1, n))
  w <- as.numeric(weights)
  if (length(w) == 1L) w <- rep(w, n)
  if (length(w) != n) stop("'weights' must have the same length as 'y'")
  if (anyNA(w) || any(!is.finite(w)) || any(w <= 0))
    stop("'weights' must be positive and finite")
  w
}

# Collapse tied values of x to their weighted mean.
#
# The Gaussian likelihood factorises as
#   sum_j w_ij (y_ij - f_i)^2 = sw_i (ybar_i - f_i)^2 + SS_i,
# so the group mean, carrying the summed weight, is a sufficient statistic for
# f.  The leftover SS_i does not involve f -- it is returned separately because
# it is exactly what identifies the observation variance, and a function that
# estimates that variance has to put it back into the likelihood.
#
# Returns the collapsed data plus the pooled within-group sum of squares and
# the number of observations it was computed from (its degrees of freedom).
ssj_collapse_ties <- function(y, x, wt) {
  if (!anyDuplicated(x))
    return(list(y = y, x = x, wt = wt, ss_within = 0, df_within = 0,
                collapsed = FALSE))
  grp <- cumsum(c(TRUE, diff(x) != 0))
  ok  <- !is.na(y)
  sw  <- as.numeric(tapply(ifelse(ok, wt, 0),     grp, sum))   # weight of the
  sy  <- as.numeric(tapply(ifelse(ok, wt*y, 0),   grp, sum))   # observed part
  ybar <- ifelse(sw > 0, sy / sw, NA_real_)
  ss   <- if (any(ok)) sum(wt[ok] * (y[ok] - ybar[grp[ok]])^2) else 0
  # a group with no observed value keeps a placeholder weight: y is NA there,
  # so the filter skips the update and the value is never used
  list(y  = ybar,
       x  = x[!duplicated(x)],
       wt = ifelse(sw > 0, sw, 1),
       ss_within = ss,
       df_within = sum(ok) - sum(sw > 0),
       collapsed = TRUE)
}

# Validate the length of a user-supplied vector of starting values.
#
# Call it AFTER ssj_collapse_ties(): the expected length is driven by the number
# of DISTINCT values of x, which with ties is not length(y), so a caller has no
# easy way of getting it right by inspection.  Without the check a wrong-length
# parinit reaches nloptr as an x0 that recycling has padded with NA, and the
# failure surfaces as the unhelpful "x0 contains NA".
#
# `layout` spells out what the vector holds, `n` and `collapsed` explain where
# the expected length comes from.
ssj_check_parinit <- function(parinit, expected, layout, n, collapsed) {
  if (is.null(parinit) || length(parinit) == expected) return(invisible(NULL))
  stop("'parinit' must have length ", expected, ", not ", length(parinit),
       ": it holds ", layout, ".  Here n = ", n,
       if (collapsed)
         paste0(", the number of distinct values of 'x' after tied values were",
                " collapsed to their weighted mean, not length(y)")
       else ", the number of observations",
       ".", call. = FALSE)
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
#' @param warp_frac share of the plotting height given to the lower panel, a
#'   number strictly between 0 and 1 (default \code{0.3}). The warp panel
#'   carries one value per interval and most of them are zero, so it needs much
#'   less room than the fitted curve; \code{0.5}, which is what equal panels
#'   amount to, wastes half the figure. It has no effect when there is no lower
#'   panel. In the \pkg{ggplot2} figure it requires \pkg{ggplot2} 4.0.0 or
#'   later, the version that introduced the \code{panel.heights} theme element;
#'   earlier versions draw panels of equal height and a non-default
#'   \code{warp_frac} then warns. Base graphics honour it at any version.
#' @param main title of the plot; the default, \code{NULL}, draws no title.
#' @param xlab,ylab axis labels; if \code{NULL} defaults are used.
#' @param col_fit,col_data colours of the fitted curve and of the observations.
#' @param use_ggplot logical: if \code{TRUE} (the default, as in
#'   \code{\link{plot.hpj}}) the figure is built with \pkg{ggplot2} and
#'   \emph{returned} rather than drawn, so it can be further modified, themed,
#'   saved with \code{ggplot2::ggsave} or combined with other plots. Set it to
#'   \code{FALSE} for base graphics, which draws directly and returns
#'   \code{NULL}. If \pkg{ggplot2} is not installed the argument is ignored
#'   with a warning and base graphics are used.
#'
#'   Note the consequence for scripts: with the default the plot is a value, so
#'   at the top level of a script or inside a loop it must be \code{print}ed to
#'   appear, exactly as any \pkg{ggplot2} figure must.
#' @param ... further arguments passed to \code{\link{ssj_breaks}}, such as
#'   \code{thresh}, \code{frac_total} or \code{frac_delta}.
#'
#' @returns With the default \code{use_ggplot = TRUE}, a \code{ggplot} object;
#'   the fitted curve and the warp profile are facets of a single plot rather
#'   than two separate figures, so that one value is returned and the two x axes
#'   stay aligned. With \code{use_ggplot = FALSE}, \code{NULL} invisibly, the
#'   function being called for its side effect of drawing.
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
#' # the default returns a ggplot object, so it can be extended
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot(fit) + ggplot2::labs(title = "Old Faithful")
#' }
#'
#' # base graphics instead
#' plot(fit, use_ggplot = FALSE)
#'
#' # give the warp panel a fifth of the height instead of the default 30%
#' plot(fit, warp_frac = 0.2)
#'
#' @export
plot.ssj <- function(x, y = NULL, prob = NULL, show_breaks = TRUE,
                     show_warp = TRUE, warp_frac = 0.3, main = NULL,
                     xlab = NULL, ylab = NULL,
                     col_fit = "firebrick", col_data = "steelblue",
                     use_ggplot = TRUE, ...) {
  obj <- x
  xx  <- obj$x
  yy  <- obj$y
  fit <- obj$smoothed_level

  if (!is.numeric(warp_frac) || length(warp_frac) != 1L ||
      !is.finite(warp_frac) || warp_frac <= 0 || warp_frac >= 1)
    stop("'warp_frac' must be a single number strictly between 0 and 1")

  if (use_ggplot && !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("package 'ggplot2' is not installed: falling back to base graphics")
    use_ggplot <- FALSE
  }

  br      <- ssj_breaks(obj, ...)
  has_warp <- (ssj_is_warp(obj) || identical(obj$method, "fix")) &&
              !is.null(if (ssj_is_warp(obj)) obj$gamma else obj$sigmas)
  panel2  <- isTRUE(show_warp) && has_warp

  # No title by default: the method labels are long, and a plot that is going
  # into a paper usually wants either no title or one the user chooses.
  if (is.null(xlab)) xlab <- "x"
  if (is.null(ylab)) ylab <- "y"

  # Confidence band, needed by both drawing paths, so computed before they
  # part.  It stays NULL when no coverage was asked for, or when the
  # estimator returns no pointwise variances (ssj_atw does not).
  band <- NULL
  if (!is.null(prob) && !is.null(obj$var_smoothed_level)) {
    z    <- abs(qnorm((1 - prob)/2))
    se   <- sqrt(pmax(obj$var_smoothed_level, 0))
    band <- list(lo = fit - z*se, up = fit + z*se)
  }

  # ---- ggplot2 branch ----------------------------------------------------
  # The two panels become facets of one plot rather than two objects, so that a
  # single value is returned and the x axes stay aligned.  Each layer is given
  # its own data frame carrying the facet key, which is how geoms that differ
  # between panels are combined without an extra package.
  if (use_ggplot) {
    key <- c("fit", "warp")
    # Short label: the panel is now about a third of the height, so a long
    # strip label crowds the axis it sits beside.
    lab_warp <- if (ssj_is_warp(obj)) "Warps" else "Additional s.d."
    lvl <- if (panel2) c(ylab_fit = ylab, w = lab_warp) else ylab
    mkp <- function(d, which) {
      d$panel <- factor(if (which == 1L) ylab else lab_warp,
                        levels = if (panel2) c(ylab, lab_warp) else ylab)
      d
    }
    gg <- ggplot2::ggplot() + ggplot2::theme_bw()

    if (!is.null(band))
      gg <- gg + ggplot2::geom_ribbon(
        data = mkp(data.frame(x = xx, lo = band$lo, up = band$up), 1L),
        ggplot2::aes(x = .data$x, ymin = .data$lo, ymax = .data$up),
        fill = col_fit, alpha = 0.20)

    if (isTRUE(show_breaks) && nrow(br) > 0L)
      gg <- gg + ggplot2::geom_vline(xintercept = br$x, linetype = "dashed",
                                     colour = "grey40")

    if (!is.null(yy))
      gg <- gg + ggplot2::geom_point(
        data = mkp(data.frame(x = xx, y = as.numeric(yy)), 1L),
        ggplot2::aes(x = .data$x, y = .data$y),
        colour = col_data, alpha = 0.35, size = 1)

    gg <- gg + ggplot2::geom_line(
      data = mkp(data.frame(x = xx, y = fit), 1L),
      ggplot2::aes(x = .data$x, y = .data$y), colour = col_fit, linewidth = 0.8)

    if (panel2) {
      g <- if (ssj_is_warp(obj)) as.numeric(obj$gamma) else as.numeric(obj$sigmas)[-1L]
      g[!is.finite(g)] <- 0
      gg <- gg + ggplot2::geom_segment(
        data = mkp(data.frame(x = xx[-1L], g = g), 2L),
        ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$g),
        colour = "grey70")
      if (isTRUE(show_breaks) && nrow(br) > 0L) {
        gg <- gg + ggplot2::geom_segment(
          data = mkp(data.frame(x = br$x, g = br$warp), 2L),
          ggplot2::aes(x = .data$x, xend = .data$x, y = 0, yend = .data$g),
          colour = col_fit, linewidth = 0.8) +
          ggplot2::geom_point(
            data = mkp(data.frame(x = br$x, g = br$warp), 2L),
            ggplot2::aes(x = .data$x, y = .data$g), colour = col_fit, size = 1.4)
      }
      gg <- gg + ggplot2::facet_grid(rows = ggplot2::vars(.data$panel),
                                     scales = "free_y", switch = "y") +
        ggplot2::theme(strip.placement = "outside",
                       strip.background = ggplot2::element_blank())
      # The warp panel carries one number per interval and needs far less room
      # than the fit.  'panel.heights' takes null units, which the gtable
      # divides in proportion, so c(1 - warp_frac, warp_frac) gives the warp
      # panel exactly that share of the plotting height.  The element arrived in
      # ggplot2 4.0.0 and older versions ignore an unknown theme element in
      # silence rather than failing, so the version is checked rather than the
      # result: an unhonoured request would otherwise look like a bug in this
      # function.  Base graphics honour 'warp_frac' at every ggplot2 version,
      # since they do not use ggplot2 at all.
      if (utils::packageVersion("ggplot2") >= "4.0.0") {
        gg <- gg + ggplot2::theme(
          panel.heights = ggplot2::unit(c(1 - warp_frac, warp_frac), "null"))
      } else if (!isTRUE(all.equal(warp_frac, 0.3))) {
        warning("'warp_frac' needs ggplot2 >= 4.0.0 to take effect in the ",
                "ggplot2 figure; the panels will be of equal height. Use ",
                "use_ggplot = FALSE for base graphics, which honour it.",
                call. = FALSE)
      }
    }
    gg <- gg + ggplot2::labs(x = xlab, y = if (panel2) NULL else ylab,
                             title = main)
    return(gg)
  }

  # ---- base graphics -----------------------------------------------------
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)

  if (panel2) {
    # layout(), not mfrow: mfrow can only make panels of equal height, and the
    # warp panel needs only 'warp_frac' of the plotting region.  Null units are
    # divided in proportion, so the two entries are the shares themselves.
    # par(op) on exit restores mfrow, which returns the device to a plain 1x1
    # layout, so no separate layout(1) is needed.
    layout(matrix(c(1L, 2L), ncol = 1L),
           heights = c(1 - warp_frac, warp_frac))
    # reserve room at the top only when there is actually a title to draw
    par(mar = c(0.5, 4.5, if (is.null(main)) 1 else 3, 1), oma = c(4, 0, 0, 0))
  }

  # ---- upper panel: data, fitted curve, declared discontinuities ----------
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
    lab <- if (ssj_is_warp(obj)) expression("Warps" ~ gamma) else "Additional s.d."

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


#' Predict from a smoothing spline with jumps
#'
#' Evaluates a fitted \code{ssj} object at arbitrary locations. For the
#' domain-warping estimators this is the cubic the model actually implies, not
#' an interpolation of the fitted values.
#'
#' The state of the warping model is the pair \eqn{(f, df/ds)} at every sampling
#' site, where \eqn{s} is the latent (warped) coordinate. Between two sites the
#' fitted function is therefore a cubic in \eqn{s}, and since \eqn{s} is a
#' linear dilation of \eqn{x} on each interval, it is a cubic in \eqn{x} as
#' well. Cubic Hermite interpolation with the stored levels and slopes
#' reproduces it exactly. On an interval carrying a large warp the same cubic is
#' traversed over a short physical distance, which is what makes the transition
#' look like a jump; a linear interpolation of the fitted values would cut that
#' corner, and a natural cubic spline through the fitted values alone would
#' smooth it away.
#'
#' Outside the range of the data the fit is extended linearly from the end knot,
#' rather than by extrapolating the end cubic.
#'
#' For an object produced by \code{\link{ssj_cssd_dp}} the call is passed to
#' \code{\link{predict_cssd_dp}}; for one produced by \code{\link{ssj_cssd}},
#' whose fit is piecewise constant by construction, a step interpolation is used.
#'
#' @param object an object of class \code{ssj}.
#' @param xnew numeric vector of locations at which to evaluate the fit; if
#'   \code{NULL} (default), the fitted values at the sampling sites are
#'   returned.
#' @param ... passed to \code{\link{predict_cssd_dp}} for CSSD objects,
#'   otherwise ignored.
#'
#' @returns Numeric vector of fitted values, in the same order as \code{xnew}.
#'
#' @seealso \code{\link{plot.ssj}}, \code{\link{ssj_breaks}}.
#'
#' @examples
#' set.seed(2026)
#' x <- sort(runif(80))
#' y <- sin(4*pi*x) + (x > 0.5) + rnorm(80, sd = 0.1)
#' # a short grid keeps the example quick; the default 21-point one is finer
#' fit <- auto_ssj_mle(y, x, ic = "ebic", ebic_xi = 0.1,
#'                     grid = seq(0, 1, length.out = 9))
#'
#' xg <- seq(0, 1, length.out = 400)
#' plot(x, y, pch = 16, cex = 0.6, col = "grey60")
#' lines(xg, predict(fit, xg), col = "firebrick", lwd = 2)
#'
#' # the cubic keeps the transition sharp where a linear interpolation of the
#' # fitted values would cut the corner
#' lines(fit$x, fit$smoothed_level, col = "steelblue", lty = 2)
#'
#' @export
predict.ssj <- function(object, xnew = NULL, ...) {
  if (is.null(xnew)) return(object$smoothed_level)
  xnew <- as.numeric(xnew)

  if (identical(object$method, "cssd"))
    return(predict_cssd_dp(object, xnew, ...))

  if (identical(object$method, "cssd_tv")) {
    # the total-variation fit is piecewise constant by construction
    return(stats::approx(object$x, object$smoothed_level, xout = xnew,
                         method = "constant", rule = 2L)$y)
  }

  xx <- object$x
  f  <- object$smoothed_level
  fp <- object$smoothed_slope
  n  <- length(xx)
  if (is.null(fp)) {
    warning("this fit carries no 'smoothed_slope': falling back to linear ",
            "interpolation. Re-fit with the current version for the cubic.")
    return(stats::approx(xx, f, xout = xnew, rule = 2L)$y)
  }
  if (n < 2L) return(rep(f[1L], length(xnew)))

  delta <- diff(xx)
  # tau = delta + gamma is the latent length of each interval; ssj_fix warps the
  # disturbance variances rather than the domain, so there tau = delta.
  gam <- if (ssj_is_warp(object)) as.numeric(object$gamma) else numeric(n - 1L)
  if (length(gam) != n - 1L) gam <- numeric(n - 1L)
  tau <- delta + gam

  i <- findInterval(xnew, xx, rightmost.closed = TRUE, all.inside = TRUE)
  u <- (xnew - xx[i]) / delta[i]          # 0..1 within the interval

  # cubic Hermite on [0, 1], with the tangents expressed over the LATENT length
  h00 <-  2*u^3 - 3*u^2 + 1
  h10 <-    u^3 - 2*u^2 + u
  h01 <- -2*u^3 + 3*u^2
  h11 <-    u^3 -   u^2
  out <- h00*f[i] + h10*tau[i]*fp[i] + h01*f[i + 1L] + h11*tau[i]*fp[i + 1L]

  # Linear extension beyond the data.  The slope is the bare latent slope, NOT
  # the physical slope just inside the end interval: the warps gamma live
  # strictly between observations, so outside the sample there is no fictitious
  # domain and the dilation is 1.  Using the inside slope would multiply it by
  # tau/delta, which on a boundary interval carrying a warp can be a factor of
  # tens and sends the extrapolation to absurd values within a short distance.
  lo <- xnew < xx[1L]
  hi <- xnew > xx[n]
  if (any(lo)) out[lo] <- f[1L] + fp[1L] * (xnew[lo] - xx[1L])
  if (any(hi)) out[hi] <- f[n]  + fp[n]  * (xnew[hi] - xx[n])
  out
}
