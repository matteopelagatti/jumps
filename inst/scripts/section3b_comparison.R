# ============================================================
# Complement to section3_comparison.R: higher noise regime.
#
# Signal: same g1 as section3:
#   g1(x) = J1(20x) + x*1[0.3,0.4](x) - x*1[0.6,1](x)
# True discontinuities at x = 0.3 (+0.3), x = 0.4 (-0.4), x = 0.6 (-0.6).
# N = 100 uniform random sites.
# Noise: sigma in {0.15, 0.25, 0.35}.
#
# What this adds relative to section3:
#
#   section3 covers sigma in {0.05, 0.10, 0.15}, where CSSD's cubic
#   smoothing spline clearly wins because the Bessel amplitude (~1.0)
#   is large relative to sigma and CSSD correctly captures the smooth
#   within-segment shape.
#
#   At sigma in {0.25, 0.35} the Bessel oscillations become harder to
#   distinguish from noise.  CSSD's cubic-spline advantage shrinks, and
#   its CV objective becomes noisier.  In this regime CSSD frequently
#   over-segments (several spurious disc detected), which inflates RMSE.
#   MLE, whose locally-constant model IS misspecified for the Bessel
#   signal, sometimes achieves lower RMSE at high sigma precisely by
#   returning fewer (or zero) detected jumps — effectively acting as a
#   regularised smoother rather than a jump detector.
#
# NOTE on ATW:
#   The auto_ssj_mle lambda estimator currently degenerates on
#   oscillatory signals (converges to lambda << 1 in most realisations).
#   A near-zero lambda feeds into auto_ssj_atw and forces its optimal
#   maxsum to zero, so ATW returns 0 detected jumps and a pure KFS
#   smooth.  This is a known limitation of the current automatic
#   parameter selection; the underlying ATW algorithm is correct and
#   works with a well-chosen fixed lambda.  Until the lambda estimation
#   is stabilised, ATW's RMSE advantage (if any) comes from
#   under-detection, not from better jump localisation.
#
# n_rep = 50 (raise to 200 for publication quality).
# ============================================================

library(jumps)

script_dir <- tryCatch(dirname(sys.frame(1L)$ofile),
                       error = function(e) "C:\\GitHub\\jumps\\inst\\scripts")

set.seed(42)

N      <- 100L
sigmas <- c(0.15, 0.25, 0.35)
n_rep  <- 50L
x_grid <- seq(0, 1, length.out = 500L)

# ---- Signal (same as section3_comparison.R) ----
g1 <- function(x) {
  besselJ(20 * x, nu = 1L) +
    x * (x >= 0.3 & x <= 0.4) -
    x * (x >= 0.6)
}
true_disc <- c(0.3, 0.4, 0.6)
true_line <- g1(x_grid)

# ---- Helpers ----
active_gamma <- function(gamma, x) {
  h <- mean(diff(sort(x)))
  gamma > pmax(max(gamma) * 0.05, 0.30 * h)
}
make_maxsum_grid <- function(y, x, n = 21L) {
  seq(0, sd(y, na.rm = TRUE) / mean(diff(x)), length.out = n)
}

# ============================================================
# Simulation
# ============================================================
cat("Starting simulation (n_rep =", n_rep, ")...\n")
results <- vector("list", length(sigmas))

for (s in seq_along(sigmas)) {
  sigma <- sigmas[s]
  cat("  sigma =", sigma, "\n")

  fits_cssd <- fits_atw <- fits_mle <- matrix(NA_real_, n_rep, length(x_grid))
  disc_cssd <- disc_atw <- disc_mle <- vector("list", n_rep)
  samp <- NULL

  for (r in seq_len(n_rep)) {
    x      <- sort(runif(N))
    y      <- g1(x) + rnorm(N, sd = sigma)
    midpts <- (x[-N] + x[-1L]) / 2L
    msg    <- make_maxsum_grid(y, x)

    # --- CSSD ---
    fit_cssd <- tryCatch(
      auto_ssj_cssd(y, x, delta = sigma),
      error = function(e) NULL
    )
    if (!is.null(fit_cssd)) {
      fits_cssd[r, ] <- predict_cssd_potts(fit_cssd, x_grid,
                                            y_orig = y, x_orig = x,
                                            delta  = sigma)
      disc_cssd[[r]] <- fit_cssd$disc_locs
    } else {
      disc_cssd[[r]] <- numeric(0L)
    }

    # --- MLE ---
    mle <- tryCatch(auto_ssj_mle(y, x, grid = msg), error = function(e) NULL)
    if (!is.null(mle)) {
      fits_mle[r, ] <- approx(mle$x, mle$smoothed_level,
                               xout = x_grid, method = "linear", rule = 2L)$y
      disc_mle[[r]] <- midpts[active_gamma(mle$gamma, x)]
    } else {
      disc_mle[[r]] <- numeric(0L)
    }

    # --- ATW (MLE-piloted lambda) ---
    # Warning: lambda is frequently near-zero for g1-type signals (see header).
    # ATW then returns 0 detected jumps and a pure smooth in most realisations.
    lam_atw <- if (!is.null(mle)) mle$pars["lambda"] else 1e-3
    atw <- tryCatch(
      auto_ssj_atw(y, x, lambda = lam_atw, grid = msg),
      error = function(e) NULL
    )
    if (!is.null(atw)) {
      fits_atw[r, ] <- approx(atw$x, atw$smoothed_level,
                               xout = x_grid, method = "linear", rule = 2L)$y
      disc_atw[[r]] <- midpts[active_gamma(atw$gamma, x)]
    } else {
      disc_atw[[r]] <- numeric(0L)
    }

    if (r == 1L)
      samp <- list(x = x, y = y, cssd = fit_cssd, atw = atw, mle = mle)

    if (r %% 10L == 0L) cat("    rep", r, "/", n_rep, "\n")
  }

  true_mat  <- matrix(true_line, nrow = n_rep, ncol = length(x_grid), byrow = TRUE)
  rmse_cssd <- sqrt(mean((fits_cssd - true_mat)^2, na.rm = TRUE))
  rmse_atw  <- sqrt(mean((fits_atw  - true_mat)^2, na.rm = TRUE))
  rmse_mle  <- sqrt(mean((fits_mle  - true_mat)^2, na.rm = TRUE))

  # Count reps where each method detected exactly the correct disc count
  n_correct_cssd <- sum(vapply(disc_cssd, length, integer(1L)) == length(true_disc), na.rm = TRUE)
  n_correct_mle  <- sum(vapply(disc_mle,  length, integer(1L)) == length(true_disc), na.rm = TRUE)
  n_correct_atw  <- sum(vapply(disc_atw,  length, integer(1L)) == length(true_disc), na.rm = TRUE)

  results[[s]] <- list(
    sigma          = sigma,
    samp           = samp,
    rmse           = c(cssd = rmse_cssd, atw = rmse_atw, mle = rmse_mle),
    n_correct      = c(cssd = n_correct_cssd, atw = n_correct_atw, mle = n_correct_mle),
    q_cssd         = apply(fits_cssd, 2L, quantile, c(.025, .975), na.rm = TRUE),
    q_atw          = apply(fits_atw,  2L, quantile, c(.025, .975), na.rm = TRUE),
    q_mle          = apply(fits_mle,  2L, quantile, c(.025, .975), na.rm = TRUE),
    disc_cssd      = disc_cssd,
    disc_atw       = disc_atw,
    disc_mle       = disc_mle
  )
  cat("    RMSE  cssd =", signif(rmse_cssd, 3),
      " atw =",  signif(rmse_atw,  3),
      " mle =",  signif(rmse_mle,  3), "\n")
  cat("    n_exact_disc cssd =", n_correct_cssd,
      " atw =", n_correct_atw, " mle =", n_correct_mle, "/", n_rep, "\n")
}

cat("Simulation complete.\n")

# ---- RMSE table ----
rmse_tab <- do.call(rbind, lapply(results, function(r) {
  c(sigma   = r$sigma,
    rmse_cssd = r$rmse["cssd"],
    rmse_atw  = r$rmse["atw"],
    rmse_mle  = r$rmse["mle"],
    ndisc_cssd = r$n_correct["cssd"],
    ndisc_atw  = r$n_correct["atw"],
    ndisc_mle  = r$n_correct["mle"])
}))
cat("\nRMSE and # reps with exactly 3 detected disc (out of", n_rep, "):\n")
print(as.data.frame(rmse_tab), digits = 4L, row.names = FALSE)

saveRDS(results, file.path(script_dir, "section3b_results.rds"))
cat("Results saved to section3b_results.rds\n")

# ============================================================
# Plotting
# ============================================================

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

draw_hist_panel <- function(disc_list, col_line, show_xlab = FALSE) {
  brk     <- seq(0, 1, by = 0.025)
  all_loc <- unlist(disc_list)
  cnt <- if (length(all_loc) > 0L)
    tabulate(findInterval(all_loc, brk, rightmost.closed = TRUE),
             nbins = length(brk) - 1L)
  else
    integer(length(brk) - 1L)
  plot(NA, xlim = c(0, 1), ylim = c(0, max(cnt, 1L)),
       xaxs = "i", yaxs = "i",
       xlab = if (show_xlab) "x" else "", ylab = "",
       main = "", axes = FALSE)
  rect(brk[-length(brk)], 0L, brk[-1L], cnt,
       col = adjustcolor(col_line, alpha.f = 0.65), border = NA)
  axis(1L, at = c(0, true_disc, 1))
  axis(2L, las = 1L, at = pretty(c(0L, max(cnt, 1L)), n = 3L))
  abline(v = true_disc, lty = 3L, col = "grey20")
}

methods   <- c("cssd",       "atw",         "mle")
meth_labs <- c("CSSD (CV)",  "ATW (MLE-λ)", "MLE")
col_lines <- c("steelblue",  "darkgreen",   "darkorange3")
sig_labs  <- paste0("σ = ", sigmas)

for (m in seq_along(methods)) {
  meth     <- methods[m]
  col_line <- col_lines[m]
  lab      <- meth_labs[m]
  fname    <- paste0("fig_g1hi_", meth, ".pdf")

  pdf(file.path(script_dir, fname), width = 10, height = 6)
  layout(matrix(seq_len(6L), nrow = 2L, ncol = 3L, byrow = TRUE),
         heights = c(3.5, 1.8))
  op <- par(mar = c(1.5, 3.5, 2.5, 0.5), mgp = c(2, 0.5, 0),
            tcl = -0.25, cex = 0.85, cex.main = 0.95)

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
    if (s == 1L) mtext("count", side = 2L, line = 2.5, cex = 0.72)
  }
  par(op)
  dev.off()
  cat("Saved", fname, "\n")
}

# Combined 6×3 figure
pdf(file.path(script_dir, "fig_g1hi_combined.pdf"), width = 10, height = 14)
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
cat("Saved fig_g1hi_combined.pdf\n")
