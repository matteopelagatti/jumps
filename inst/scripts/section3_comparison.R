# ============================================================
# Replication of the simulated experiment from Section 3 of
#   Storath & Weinmann (2024), "Smoothing splines for
#   discontinuous signals", JCGS 33(2):651-664.
#
# Three fully-automatic methods are compared, each producing
# its own 2 x 3 PDF figure (3 sigma levels, curve + histogram):
#
#   CSSD-CV  auto_ssj_cssd: 5-fold CV selects (p, gamma) jointly;
#            implements the exact weak-rod model from Storath & Weinmann (2024)
#   ATW      auto_ssj_atw with lambda from MLE pilot, EBIC selects maxsum
#   MLE      auto_ssj_mle: MLE selects (sigma, lambda), EBIC selects maxsum
#
# Signal:  g1(x) = J1(20x) + x*1[0.3,0.4](x) - x*1[0.6,1](x)
#          (J1 = Bessel function of the first kind, order 1)
# True discontinuities at x = 0.3 (+0.3), x = 0.4 (-0.4), x = 0.6 (-0.6).
# N = 100 uniform random sites in [0, 1].  Noise: sigma in {0.05, 0.10, 0.15}.
#
# NOTE: ssj_cssd requires the package to be rebuilt after the unweighted-D
# fix in src/cssd.cpp (remove the 1/h weighting that caused ADMM to stall).
# Run devtools::install("path/to/jumps") after rebuilding.
# ============================================================

library(jumps)

script_dir <- tryCatch(dirname(sys.frame(1L)$ofile), error = function(e) "C:\\GitHub\\jumps\\inst\\scripts")

set.seed(2024)

# ---- Adjustable parameters ----
N      <- 100L
sigmas <- c(0.05, 0.10, 0.15)
n_rep  <- 100L     # raise to 500 for publication quality; 100 takes ~20-30 min
x_grid <- seq(0, 1, length.out = 500L)

# ---- Signal ----
# g1(x) = J1(20x) + x*1[0.3,0.4](x) - x*1[0.6,1](x)
g1 <- function(x) {
  besselJ(20 * x, nu = 1L) +
    x * (x >= 0.3 & x <= 0.4) -
    x * (x >= 0.6)
}
true_disc <- c(0.3, 0.4, 0.6)
true_line <- g1(x_grid)

# ====================================================================
# Helper: 5-fold CV lambda selection for ssj_cssd (TV denoising)
# ====================================================================
# The TV penalty with unweighted D (entries +/-1) uses lambda directly as
# the soft threshold on absolute first differences.  The smooth Bessel signal
# contributes differences of order max|g1'|*h ~ 0.2 (for h~0.01), which
# overlaps the jump amplitudes 0.3-0.6.  The fixed-formula "2*sigma" rule
# is therefore too tight; k-fold CV learns the right threshold from data.
#
# Lambda grid: spans from just below the noise floor (sigma*sqrt(2)~0.07)
# to well above the largest jump (0.6), on a log scale.
lambda_grid_tv <- exp(seq(log(0.02), log(1.5), length.out = 25L))

cssd_cv_lambda <- function(y, x,
                            lambda_grid = lambda_grid_tv,
                            k_folds     = 5L, ...) {
  n       <- length(y)
  fold_id <- sample(rep(seq_len(k_folds), length.out = n))
  cv_err  <- vapply(lambda_grid, function(lam) {
    err <- 0
    for (fold in seq_len(k_folds)) {
      tr     <- which(fold_id != fold)
      te     <- which(fold_id == fold)
      ord_tr <- order(x[tr])
      fit  <- ssj_cssd(y[tr][ord_tr], x[tr][ord_tr], lambda = lam, ...)
      pred <- approx(fit$x, fit$smoothed_level,
                     xout = x[te], method = "constant", rule = 2L)$y
      err  <- err + sum((y[te] - pred)^2, na.rm = TRUE)
    }
    err
  }, numeric(1L))
  lambda_grid[which.min(cv_err)]
}

# ====================================================================
# Helper: maxsum grid for auto methods
# ====================================================================
# gamma_t (the virtual time added at gap t) is in x-units.  In the KFS/ATW
# model with lambda~1e-4 and h~0.01, creating a visible jump requires
# gamma_t ~ several times delta.  sd(y)/mean(diff(x)) ~ 40 is a generous
# upper bound; early stopping in auto_ssj_* means at most 7-10 of these
# values are actually evaluated per replication.
make_maxsum_grid <- function(y, x, n = 21L) {
  seq(0, sd(y, na.rm = TRUE) / mean(diff(x)), length.out = n)
}

# ====================================================================
# Jump-threshold helper
# ====================================================================
# The L2 budget constraint (sum(gamma) <= maxsum) does NOT induce sparsity:
# the optimizer can distribute budget across many small gamma_t values.
# Two complementary thresholds are used:
#   (a) relative: must exceed 5% of the largest gamma_t (catches clusters
#       of similarly-sized spurious values when max is itself real);
#   (b) absolute: must exceed 30% of the mean inter-observation spacing h
#       (a genuine time-warp requires gamma_t to be meaningfully large
#       relative to the interval being stretched).
# A gamma_t must clear BOTH thresholds to be flagged as an active jump.
active_gamma <- function(gamma, x) {
  h   <- mean(diff(sort(x)))
  rel <- max(gamma) * 0.05          # (a) 5 % of the largest gamma
  abs <- 0.30 * h                   # (b) 30 % of mean spacing
  gamma > pmax(rel, abs)
}

# ====================================================================
# Simulation
# ====================================================================

cat("Starting simulation (n_rep =", n_rep, ")...\n")

results <- vector("list", length(sigmas))

for (s in seq_along(sigmas)) {
  sigma <- sigmas[s]
  cat("  sigma =", sigma, "\n")

  fits_tv  <- fits_atw <- fits_mle <- matrix(NA_real_, n_rep, length(x_grid))
  disc_tv  <- disc_atw <- disc_mle <- vector("list", n_rep)
  samp <- NULL

  for (r in seq_len(n_rep)) {
    x      <- sort(runif(N))
    y      <- g1(x) + rnorm(N, sd = sigma)
    midpts <- (x[-N] + x[-1L]) / 2L
    msg    <- make_maxsum_grid(y, x)

    # --- CSSD (true cubic smoothing spline with Potts penalty, auto CV) ---
    fit_cssd <- tryCatch(
      auto_ssj_cssd(y, x, delta = sigma),
      error = function(e) NULL
    )
    if (!is.null(fit_cssd)) {
      fits_tv[r, ] <- predict_cssd_potts(fit_cssd, x_grid,
                                          y_orig = y, x_orig = x,
                                          delta = sigma)
      disc_tv[[r]] <- fit_cssd$disc_locs
    } else {
      disc_tv[[r]] <- numeric(0L)
    }

    # --- MLE (auto_ssj_mle): run first so lambda can be reused by ATW ---
    mle <- tryCatch(
      auto_ssj_mle(y, x),
      error = function(e) NULL
    )
    if (!is.null(mle)) {
      fits_mle[r, ] <- approx(mle$x, mle$smoothed_level,
                               xout = x_grid, method = "linear", rule = 2L)$y
      disc_mle[[r]] <- midpts[active_gamma(mle$gamma, x)]
    }

    # --- ATW (MLE-piloted lambda, EBIC selects maxsum) ---
    # GCV is unreliable for jump signals: a flexible spline (small lambda)
    # can explain discontinuities through curvature, giving near-zero GCV
    # score at the wrong lambda.  Using the MLE-estimated lambda avoids this
    # bias; the maxsum budget is still selected by EBIC via auto_ssj_atw.
    lam_atw <- if (!is.null(mle)) mle$pars["lambda"] else 1e-3
    atw <- tryCatch(
      auto_ssj_atw(y, x, lambda = lam_atw),
      error = function(e) NULL
    )
    if (!is.null(atw)) {
      fits_atw[r, ] <- approx(atw$x, atw$smoothed_level,
                               xout = x_grid, method = "linear", rule = 2L)$y
      disc_atw[[r]] <- midpts[active_gamma(atw$gamma, x)]
    }

    if (r == 1L)
      samp <- list(x = x, y = y, tv = fit_cssd, atw = atw, mle = mle)

    if (r %% 20L == 0L) cat("    rep", r, "/", n_rep, "\n")
  }

  true_mat <- matrix(true_line, nrow = n_rep, ncol = length(x_grid), byrow = TRUE)
  rmse_tv  <- sqrt(mean((fits_tv  - true_mat)^2, na.rm = TRUE))
  rmse_atw <- sqrt(mean((fits_atw - true_mat)^2, na.rm = TRUE))
  rmse_mle <- sqrt(mean((fits_mle - true_mat)^2, na.rm = TRUE))

  results[[s]] <- list(
    sigma = sigma,
    samp  = samp,
    rmse  = c(cssd = rmse_tv, atw = rmse_atw, mle = rmse_mle),
    q_tv  = apply(fits_tv,  2L, quantile, c(.025, .975), na.rm = TRUE),
    q_atw = apply(fits_atw, 2L, quantile, c(.025, .975), na.rm = TRUE),
    q_mle = apply(fits_mle, 2L, quantile, c(.025, .975), na.rm = TRUE),
    disc_tv  = disc_tv,
    disc_atw = disc_atw,
    disc_mle = disc_mle
  )
  cat("    CSSD done (sigma =", sigma, ")\n")
}

cat("Simulation complete.\n")

# ---- RMSE table ----
rmse_tab <- do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$rmse)))
cat("\nRMSE (signal estimate vs. true g1, averaged over reps and x_grid):\n")
print(as.data.frame(rmse_tab), digits = 4L, row.names = FALSE)

saveRDS(results, file.path(script_dir, "section3_results.rds"))
cat("Results saved to section3_results.rds\n")

# ====================================================================
# Plotting utilities
# ====================================================================

# Fitted-curve panel: data + true signal + MC band + sample fit + jump marks
draw_fit_panel <- function(res, method, col_line,
                            main = "", ylab = "", show_xlab = FALSE) {
  smp      <- res$samp
  fit      <- smp[[method]]
  q        <- res[[paste0("q_", method)]]
  col_band <- adjustcolor(col_line, alpha.f = 0.20)
  ylim     <- range(c(smp$y, q, true_line), na.rm = TRUE)
  ylim     <- ylim + c(-1, 1) * 0.06 * diff(ylim)

  plot(smp$x, smp$y, pch = 1L, col = "grey55", cex = 0.5,
       xlim = c(0, 1), ylim = ylim,
       xlab = if (show_xlab) "x" else "",
       ylab = ylab, main = main, xaxs = "i", las = 1L)
  lines(x_grid, true_line, lty = 2L, col = "grey40", lwd = 0.9)
  polygon(c(x_grid, rev(x_grid)), c(q[1L, ], rev(q[2L, ])),
          col = col_band, border = NA)
  if (!is.null(fit))
    lines(fit$x, fit$smoothed_level, col = col_line, lwd = 1.8)
  abline(v = true_disc, lty = 3L, col = "grey20")
}

# Histogram panel: frequency of detected jump locations across reps
draw_hist_panel <- function(disc_list, col_line, show_xlab = FALSE) {
  brk     <- seq(0, 1, by = 0.025)
  all_loc <- unlist(disc_list)
  cnt <- if (length(all_loc) > 0L) {
    tabulate(findInterval(all_loc, brk, rightmost.closed = TRUE),
             nbins = length(brk) - 1L)
  } else {
    integer(length(brk) - 1L)
  }
  plot(NA, xlim = c(0, 1), ylim = c(0, max(cnt, 1L)),
       xaxs = "i", yaxs = "i",
       xlab = if (show_xlab) "x" else "", ylab = "",
       main = "", axes = FALSE)
  rect(brk[-length(brk)], 0L, brk[-1L], cnt,
       col = adjustcolor(col_line, alpha.f = 0.65), border = NA)
  axis(1L, at = c(0, 0.3, 0.4, 0.6, 1))
  axis(2L, las = 1L, at = pretty(c(0L, max(cnt, 1L)), n = 3L))
  abline(v = true_disc, lty = 3L, col = "grey20")
}

# ====================================================================
# Figures: one 2 x 3 PDF per method
# ====================================================================

methods   <- c("tv",           "atw",              "mle")
meth_labs <- c("CSSD (CV)",    "ATW (MLE-λ)",      "MLE")
col_lines <- c("steelblue",    "darkgreen",         "darkorange3")
sig_labs  <- paste0("σ = ", sigmas)

for (m in seq_along(methods)) {
  meth     <- methods[m]
  col_line <- col_lines[m]
  lab      <- meth_labs[m]
  fname    <- paste0("fig_", meth, ".pdf")

  pdf(file.path(script_dir, fname), width = 10, height = 6)
  layout(matrix(seq_len(6L), nrow = 2L, ncol = 3L, byrow = TRUE),
         heights = c(3.5, 1.8))
  op <- par(mar = c(1.5, 3.5, 2.5, 0.5), mgp = c(2, 0.5, 0),
            tcl = -0.25, cex = 0.85, cex.main = 0.95)

  # Row 1: sample fit + MC confidence band
  for (s in seq_along(sigmas)) {
    draw_fit_panel(results[[s]], method = meth, col_line = col_line,
                   main  = paste0(lab, "  (", sig_labs[s], ")"),
                   ylab  = if (s == 1L) "y" else "")
  }

  # Row 2: histogram of detected jump locations
  par(mar = c(3, 3.5, 0.3, 0.5))
  for (s in seq_along(sigmas)) {
    draw_hist_panel(results[[s]][[paste0("disc_", meth)]],
                    col_line  = col_line,
                    show_xlab = TRUE)
    if (s == 1L) mtext("count", side = 2L, line = 2.5, cex = 0.72)
  }
  par(op)
  dev.off()
  cat("Saved", fname, "\n")
}

# ====================================================================
# Optional: combined 6 x 3 figure showing all three methods together
# ====================================================================
# Layout: rows alternate curve / histogram for each method (6 rows, 3 cols)

pdf(file.path(script_dir, "fig_combined.pdf"), width = 10, height = 14)
mat <- matrix(seq_len(18L), nrow = 6L, ncol = 3L, byrow = TRUE)
layout(mat, heights = rep(c(3.5, 1.8), 3L))
op <- par(mar = c(1.5, 3.5, 2.5, 0.5), mgp = c(2, 0.5, 0),
          tcl = -0.25, cex = 0.80, cex.main = 0.93)

for (m in seq_along(methods)) {
  meth     <- methods[m]
  col_line <- col_lines[m]
  lab      <- meth_labs[m]

  for (s in seq_along(sigmas)) {
    draw_fit_panel(results[[s]], method = meth, col_line = col_line,
                   main = paste0(lab, "  (", sig_labs[s], ")"),
                   ylab = if (s == 1L) "y" else "")
  }

  par(mar = c(3, 3.5, 0.3, 0.5))
  for (s in seq_along(sigmas)) {
    draw_hist_panel(results[[s]][[paste0("disc_", meth)]],
                    col_line  = col_line,
                    show_xlab = TRUE)
    if (s == 1L) mtext("count", side = 2L, line = 2.5, cex = 0.70)
  }
  par(mar = c(1.5, 3.5, 2.5, 0.5))
}

par(op)
dev.off()
cat("Saved fig_combined.pdf\n")
