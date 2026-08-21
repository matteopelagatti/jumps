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
#   ATW-CV   auto_ssj_atw_cv: 5-fold CV selects lambda, EBIC selects maxsum.
#            Self-contained, unlike the earlier MLE-piloted variant: the fast
#            algorithm no longer needs the estimator it is meant to
#            approximate in order to run, and it is tuned by cross-validation
#            on the same footing as CSSD.
#   MLE      auto_ssj_mle: MLE selects (sigma, lambda), EBIC selects maxsum
#
# Timing: each method is timed separately inside every replication, in CPU
# seconds (user + system) rather than wall clock, so the figures do not
# depend on how many replications are running concurrently.
#
# Signal:  g1(x) = J1(20x) + x*1[0.3,0.4](x) - x*1[0.6,1](x)
#          (J1 = Bessel function of the first kind, order 1)
# True discontinuities at x = 0.3 (+0.3), x = 0.4 (-0.4), x = 0.6 (-0.6).
# N = 100 uniform random sites in [0, 1].  Noise: sigma in {0.05, 0.10, 0.15}.
#
# Evaluation rule: every estimator is scored on x_grid through its OWN exact
# predictor - predict() on the fitted object - rather than by interpolating its
# fitted values.  For the domain-warping estimators that is the cubic the state
# (f, df/ds) implies; for CSSD the piecewise cubic of the weak-rod fit; for the
# smoothing spline predict.smooth.spline().  Scoring one method with its exact
# curve and another with a linear interpolation would favour the latter across
# the jumps.  The grid is clamped to the observed range of x before evaluation,
# so that no estimator is scored on extrapolation (see run_rep below).
#
# Common random numbers: replication r uses the same design and the same
# standardised noise at all three noise levels.  Each replication seeds its own
# RNG from seed_base + r, so the results do not depend on how the work is split
# across cores: a parallel run reproduces a serial one exactly.
# ============================================================

library(jumps)

script_dir <- tryCatch(dirname(sys.frame(1L)$ofile), error = function(e) "C:\\GitHub\\jumps\\inst\\scripts")

# One base seed.  Replication r is generated from seed_base + r at every noise
# level (common random numbers), so the seed is also set inside the loop.
seed_base <- 2024L
set.seed(seed_base)

# ---- Adjustable parameters ----
N      <- 100L
sigmas <- c(0.05, 0.10, 0.15)
n_rep  <- 1000L
x_grid <- seq(0, 1, length.out = 500L)

# Cores used for the replication loop.  Set to 1 to run serially (useful when
# debugging: errors inside a worker are harder to read).  Measured on a 4-core
# (8-thread) 2026 laptop at N = 100: about 0.56 s per replication serially, so
# the full 3 x 1000 design takes roughly 28 min on one core and 11 min on seven.
# The speedup is well below the core count because only 4 cores are physical.
n_cores <- max(1L, parallel::detectCores() - 1L)

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
# Jump-threshold helper
# ====================================================================
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

# One replication: everything that depends on r and sigma, and nothing else.
# It re-seeds from seed_base + r, so it is self-contained and its result does
# not depend on the order in which replications are executed -- which is what
# makes the parallel and the serial runs agree exactly.
run_rep <- function(r, sigma) {
  set.seed(seed_base + r)
  x      <- sort(runif(N))
  z      <- rnorm(N)
  y      <- g1(x) + sigma * z
  midpts <- (x[-N] + x[-1L]) / 2L

  # The evaluation grid spans [0, 1] but the design does not: with N uniform
  # draws roughly 1% of the grid falls outside [min(x), max(x)] at each end.
  # Scoring an estimator there measures extrapolation, not signal recovery, and
  # it is wildly unstable: a discontinuity near the boundary can leave a
  # terminal segment of two points whose end derivative is essentially
  # unconstrained, and a linear extension of it produces errors of several
  # units over a gap of a few hundredths.  Every estimator is therefore
  # evaluated at the grid CLAMPED to the observed range, which is the
  # rule = 2 convention the earlier version of this script used for the
  # time-warping estimators, now applied uniformly to all four.
  x_eval <- pmin(pmax(x_grid, min(x)), max(x))

  out <- list(cssd = rep(NA_real_, length(x_grid)),
              mle  = rep(NA_real_, length(x_grid)),
              atw  = rep(NA_real_, length(x_grid)),
              ss   = rep(NA_real_, length(x_grid)),
              disc_cssd = numeric(0L), disc_mle = numeric(0L),
              disc_atw  = numeric(0L),
              fail = c(cssd = 0L, atw = 0L, mle = 0L, ss = 0L),
              time = c(cssd = NA_real_, atw = NA_real_,
                       mle  = NA_real_, ss  = NA_real_),
              samp = NULL)

  # CPU seconds, so a replication timed while six others share the machine is
  # comparable with one timed alone.  The argument is forced inside, after the
  # clock is read, so the whole fit is inside the interval.
  timed <- function(expr) {
    t0  <- proc.time()
    val <- expr
    dt  <- proc.time() - t0
    list(val = val, secs = unname(dt[["user.self"]] + dt[["sys.self"]]))
  }

  # --- CSSD (exact weak-rod model, 5-fold CV over (p, gamma)) ---
  tc <- timed(tryCatch(auto_ssj_cssd(y, x, delta = sigma), error = function(e) NULL))
  fit_cssd <- tc$val; out$time["cssd"] <- tc$secs
  if (!is.null(fit_cssd)) {
    out$cssd      <- predict(fit_cssd, x_eval)
    out$disc_cssd <- fit_cssd$disc_locs
  } else out$fail["cssd"] <- 1L

  # --- MLE ---
  tc <- timed(tryCatch(auto_ssj_mle(y, x, ic = "ebic", ebic_xi = 0.3),
                       error = function(e) NULL))
  mle <- tc$val; out$time["mle"] <- tc$secs
  if (!is.null(mle)) {
    out$mle      <- predict(mle, x_eval)
    out$disc_mle <- midpts[active_gamma(mle$gamma, x)]
  } else out$fail["mle"] <- 1L

  # --- ATW (lambda by 5-fold CV, EBIC selects maxsum) ---
  # auto_ssj_atw_cv() cross-validates over (lambda, maxsum), keeps only lambda,
  # and hands it to auto_ssj_atw() for the budget.  The budget is deliberately
  # NOT taken from the cross-validation: prediction error is close to
  # indifferent to a few spurious warps, so it over-declares discontinuities.
  # The boundary warning is suppressed here because a grid endpoint is reported
  # per replication and would flood the console; the selected lambda is stored
  # instead, so the boundary rate can be checked afterwards.
  tc <- timed(tryCatch(suppressWarnings(
          auto_ssj_atw_cv(y, x, ic = "ebic", ebic_xi = 0.3)),
          error = function(e) NULL))
  atw <- tc$val; out$time["atw"] <- tc$secs
  if (!is.null(atw)) {
    out$atw      <- predict(atw, x_eval)
    out$disc_atw <- midpts[active_gamma(atw$gamma, x)]
  } else out$fail["atw"] <- 1L

  # --- Smoothing spline (GCV-selected) ---
  tc <- timed(tryCatch(smooth.spline(x, y), error = function(e) NULL))
  fit_ss <- tc$val; out$time["ss"] <- tc$secs
  if (!is.null(fit_ss)) out$ss <- predict(fit_ss, x_eval)$y
  else out$fail["ss"] <- 1L

  # keep the first replication's fitted objects for the illustrative panels
  if (r == 1L)
    out$samp <- list(x = x, y = y, cssd = fit_cssd, atw = atw, mle = mle, ss = fit_ss)
  out
}

# ---- worker pool ---------------------------------------------------------
# Replications are independent, so the r-loop parallelises directly.  Each
# replication seeds its own RNG, so no cluster RNG stream is needed and the
# results are identical to a serial run.
cl <- NULL
if (n_cores > 1L) {
  cl <- parallel::makeCluster(n_cores)
  parallel::clusterEvalQ(cl, library(jumps))
  parallel::clusterExport(cl,
    c("g1", "N", "x_grid", "seed_base", "active_gamma"), envir = environment())
}
# NB: no on.exit() here.  At the top level of a script it is a no-op under
# Rscript, but under source() it attaches to the frame of the current top-level
# expression and would stop the cluster immediately.  The cluster is shut down
# in the finally block below instead.

cat("Starting simulation (n_rep =", n_rep, "on",
    if (is.null(cl)) 1L else n_cores, "core(s))...\n")

results <- vector("list", length(sigmas))

tryCatch({
for (s in seq_along(sigmas)) {
  sigma <- sigmas[s]
  cat("  sigma =", sigma, "... ")
  t0 <- proc.time()

  reps <- if (is.null(cl)) {
    lapply(seq_len(n_rep), run_rep, sigma = sigma)
  } else {
    parallel::parLapplyLB(cl, seq_len(n_rep), run_rep, sigma = sigma)
  }

  grab <- function(nm) do.call(rbind, lapply(reps, `[[`, nm))
  fits_cssd <- grab("cssd"); fits_mle <- grab("mle")
  fits_atw  <- grab("atw");  fits_ss  <- grab("ss")
  disc_cssd <- lapply(reps, `[[`, "disc_cssd")
  disc_mle  <- lapply(reps, `[[`, "disc_mle")
  disc_atw  <- lapply(reps, `[[`, "disc_atw")
  fail      <- Reduce(`+`, lapply(reps, `[[`, "fail"))
  samp      <- reps[[1L]]$samp
  # CPU seconds per fit, per method.  The whole n_rep x 4 matrix is kept, not
  # just its summaries: the distribution is right-skewed (a handful of fits take
  # several times the median), so a mean alone is misleading and a box plot or
  # any other quantile can be drawn afterwards from the stored object.
  tmat      <- do.call(rbind, lapply(reps, `[[`, "time"))
  t_mean    <- colMeans(tmat, na.rm = TRUE)
  t_med     <- apply(tmat, 2L, median, na.rm = TRUE)
  t_qs      <- apply(tmat, 2L, quantile, c(.25, .75), na.rm = TRUE)

  cat(sprintf("%.1f min\n", (proc.time() - t0)[3]/60))

  true_mat <- matrix(true_line, nrow = n_rep, ncol = length(x_grid), byrow = TRUE)

  # Per-replication RMSE, so a Monte Carlo standard error can be attached to the
  # pooled figure.  Failed fits leave a row of NAs: they are dropped here, but
  # their number is reported rather than passing silently.
  rmse_rep <- function(fits) sqrt(rowMeans((fits - true_mat)^2))
  rr <- list(cssd = rmse_rep(fits_cssd), atw = rmse_rep(fits_atw),
             mle  = rmse_rep(fits_mle),  ss  = rmse_rep(fits_ss))
  mc_se <- vapply(rr, function(v) {
    v <- v[is.finite(v)]
    if (length(v) < 2L) NA_real_ else sd(v)/sqrt(length(v))
  }, numeric(1))
  # Keep every per-replication RMSE, as for the timings.  The distribution has
  # a long right tail -- CSSD in particular diverges on a couple of draws in a
  # thousand, when a discontinuity isolates two nearly coincident design points
  # into a segment of their own -- so the pooled figure and the median tell
  # different stories, and both belong in the record.
  rmse_all  <- do.call(cbind, rr)
  rmse_med  <- vapply(rr, function(v) median(v[is.finite(v)]), numeric(1))
  if (any(fail > 0L))
    cat("    failed fits:",
        paste(names(fail)[fail > 0L], fail[fail > 0L], sep = "=", collapse = ", "),
        "\n")

  rmse_cssd <- sqrt(mean((fits_cssd - true_mat)^2, na.rm = TRUE))
  rmse_atw  <- sqrt(mean((fits_atw  - true_mat)^2, na.rm = TRUE))
  rmse_mle  <- sqrt(mean((fits_mle  - true_mat)^2, na.rm = TRUE))
  rmse_ss   <- sqrt(mean((fits_ss   - true_mat)^2, na.rm = TRUE))

  results[[s]] <- list(
    sigma     = sigma,
    samp      = samp,
    rmse      = c(cssd = rmse_cssd, atw = rmse_atw, mle = rmse_mle, ss = rmse_ss),
    rmse_se   = mc_se,
    rmse_med  = rmse_med,
    rmse_all  = rmse_all,      # n_rep x 4: every per-replication RMSE kept
    n_fail    = fail,
    n_disc    = c(cssd = mean(lengths(disc_cssd)), atw = mean(lengths(disc_atw)),
                  mle  = mean(lengths(disc_mle))),
    time_mean = t_mean,
    time_med  = t_med,
    time_iqr  = t_qs,
    time_all  = tmat,          # n_rep x 4, CPU seconds: every measurement kept
    mean_cssd = colMeans(fits_cssd, na.rm = TRUE),
    mean_atw  = colMeans(fits_atw,  na.rm = TRUE),
    mean_mle  = colMeans(fits_mle,  na.rm = TRUE),
    mean_ss   = colMeans(fits_ss,   na.rm = TRUE),
    mse_cssd  = colMeans((fits_cssd - true_mat)^2, na.rm = TRUE),
    mse_atw   = colMeans((fits_atw  - true_mat)^2, na.rm = TRUE),
    mse_mle   = colMeans((fits_mle  - true_mat)^2, na.rm = TRUE),
    mse_ss    = colMeans((fits_ss   - true_mat)^2, na.rm = TRUE),
    q_cssd    = apply(fits_cssd, 2L, quantile, c(.025, .975), na.rm = TRUE),
    q_atw     = apply(fits_atw,  2L, quantile, c(.025, .975), na.rm = TRUE),
    q_mle     = apply(fits_mle,  2L, quantile, c(.025, .975), na.rm = TRUE),
    q_ss      = apply(fits_ss,   2L, quantile, c(.025, .975), na.rm = TRUE),
    disc_cssd = disc_cssd,
    disc_atw  = disc_atw,
    disc_mle  = disc_mle
  )
}
}, finally = if (!is.null(cl)) { parallel::stopCluster(cl); cl <- NULL })

cat("Simulation complete.\n")

# ---- RMSE table ----------------------------------------------------------
# CSSD is reported twice.  On a few replications its fit diverges BETWEEN the
# design points: when a declared discontinuity isolates two nearly coincident
# design points into a segment of their own, the smoothing spline on that
# segment is just the straight line through the two of them, and the linear
# extension of Storath and Weinmann's linext_pp carries that slope -- of order
# 1000 here -- across the segment's whole domain.  The fitted values AT the
# design points stay sane; only the evaluation between them explodes.  This is
# faithful to the reference implementation, not an artefact of ours, but four
# replications out of 3000 then move the pooled RMSE by up to 41%, so both
# figures are given.
#
# A replication counts as divergent when its RMSE exceeds DIVERGENT_MULT times
# the median RMSE of the same method at the same noise level.  The rule is
# scale free, and the gap it has to straddle here is wide: the flagged
# replications sit at 8-23 times the median, the largest unflagged one at
# about 3.7, so the count is insensitive to the multiplier.
DIVERGENT_MULT <- 5

# The table is built from `rmse_all`, which the per-sigma loop above fills in.
# Guard against reaching this point with a `results` left over from an older
# version of the script that did not store it: NULL[, "cssd"] is NULL rather
# than an error, so without this check the table would quietly print NaN
# instead of failing.  In a run from the top the condition never fires.
if (!all(vapply(results, function(r) is.matrix(r$rmse_all), TRUE)))
  stop("`results` carries no 'rmse_all': it predates this table. ",
       "Re-run the script from the top.", call. = FALSE)

divergent <- function(v, mult = DIVERGENT_MULT) {
  v  <- as.numeric(v)
  ok <- is.finite(v)
  ok & v > mult * median(v[ok])
}

# Pooled RMSE over replications and grid points; identical to
# sqrt(mean((fits - truth)^2)) because every replication shares the same grid.
# An empty selection returns NA rather than NaN, so a column that ends up with
# nothing in it reads as missing instead of as a number.
pooled <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) NA_real_ else sqrt(mean(v^2))
}

rmse_tab <- do.call(rbind, lapply(results, function(r) {
  v   <- r$rmse_all
  bad <- divergent(v[, "cssd"])
  c(sigma                  = r$sigma,
    "CSSD (all)"           = pooled(v[, "cssd"]),
    "CSSD (excl. diverg.)" = pooled(v[!bad, "cssd"]),
    "SSJ-ATW"              = pooled(v[, "atw"]),
    "SSJ-MLE"              = pooled(v[, "mle"]),
    "Smooth. spline"       = pooled(v[, "ss"]),
    "n excl."              = sum(bad))
}))

cat("\nRMSE (signal estimate vs. true g1, over replications and x_grid).\n")
cat("CSSD is shown with and without the replications on which its fit\n")
cat("diverges between design points; 'n excl.' counts those.\n")
print(as.data.frame(rmse_tab), digits = 4L, row.names = FALSE)

# ---- LaTeX version of the RMSE table --------------------------------------
# Emitted from the same rmse_tab that was just printed, so the file and the
# console cannot disagree, and the paper cannot drift from the simulation.
# Values are scaled by 1000 and given to one decimal, the form used in the
# paper; the caption says so, because g1 ranges over roughly [-1.1, 0.6] and a
# reader who missed the scaling would place the errors two orders of magnitude
# out.  The replication and grid counts are interpolated rather than written
# in, so a change to n_rep cannot leave a stale number in the caption.
#
# The count of excluded replications sits inside the CSSD group rather than at
# the end of the row: it qualifies those two columns and nothing else, and a
# lone trailing column of mostly zeros invites the reader to take it as a
# property of the design instead.
tex_file <- file.path(script_dir, "tab_rmse.tex")

tex_body <- vapply(seq_len(nrow(rmse_tab)), function(i) {
  z <- rmse_tab[i, ]
  sprintf("%.2f & %5.1f & %5.1f & %d & %5.1f & %5.1f & %5.1f \\\\",
          z[["sigma"]],
          1000 * z[["CSSD (all)"]],
          1000 * z[["CSSD (excl. diverg.)"]],
          as.integer(z[["n excl."]]),
          1000 * z[["SSJ-ATW"]],
          1000 * z[["SSJ-MLE"]],
          1000 * z[["Smooth. spline"]])
}, character(1))

writeLines(c(
  "% Generated by inst/scripts/section3_comparison.R -- do not edit by hand.",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Root mean squared error of the fitted curve against the true",
  sprintf("signal $g_1$, pooled over %d replications and the %d-point", n_rep,
          length(x_grid)),
  "evaluation grid, at three noise levels. All values are multiplied by",
  "$10^3$. CSSD is reported twice: over all replications, and excluding the",
  "$n_{\\mathrm{ex}}$ replications on which its fit diverges between the",
  "design points.}",
  "\\label{tab:rmse}",
  "\\begin{tabular}{lccc ccc}",
  "\\toprule",
  " & \\multicolumn{3}{c}{CSSD} & \\multicolumn{3}{c}{} \\\\",
  "\\cmidrule(lr){2-4}",
  paste("$\\sigma$ & all & excl.\\ divergent & $n_{\\mathrm{ex}}$ &",
        "SSJ--ATW & SSJ--MLE & Smoothing spline \\\\"),
  "\\midrule",
  tex_body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), tex_file)
cat("Saved", basename(tex_file), "\n")


# ---- Detection table ------------------------------------------------------
# Built from the same `results` as the RMSE table, and written to disk in the
# same pass, so the two tables in the paper cannot come from different runs.
# (They did once: an earlier detection table survived the ebic_xi calibration
# that changed how many discontinuities the estimators declare, and reported
# SSJ-MLE detecting far less often than SSJ-ATW -- an ordering that no longer
# holds.)
#
# A break counts as detected in a replication when at least one discontinuity
# is flagged within DETECT_TOL of it; the rate is the fraction of replications,
# not the number of flags, so a break declared twice in one replication counts
# once.
#
# The denominator for the false positive rate is NOT the evaluation grid.
# Discontinuities can only be declared at the N-1 midpoints between consecutive
# design points, so that is the candidate set; of those, the ones outside all
# three windows are the opportunities to raise a false alarm.  Their number
# varies a little with the design, so it is averaged over replications and
# reported in the caption.  The raw count per replication is given as well,
# because it needs no denominator at all.
DETECT_TOL <- 0.02

in_window <- function(v) {
  if (!length(v)) return(logical(0L))
  Reduce(`|`, lapply(true_disc, function(b) abs(v - b) <= DETECT_TOL))
}

# Candidate midpoints outside the windows, regenerated from the same seeds as
# the simulation itself, so the denominator matches the designs actually used.
n_out_mean <- mean(vapply(seq_len(n_rep), function(r) {
  set.seed(seed_base + r)
  x <- sort(runif(N))
  sum(!in_window((x[-N] + x[-1L]) / 2))
}, numeric(1L)))

det_row <- function(d) {
  det <- vapply(true_disc,
                function(b) mean(vapply(d, function(v) any(abs(v - b) <= DETECT_TOL),
                                        logical(1L))),
                numeric(1L))
  fp <- mean(vapply(d, function(v) sum(!in_window(v)), numeric(1L)))
  c(det, fp = fp, fp_rate = fp / n_out_mean)
}

meth_lab <- c(disc_cssd = "CSSD", disc_atw = "SSJ--ATW", disc_mle = "SSJ--MLE")

det_tab <- do.call(rbind, lapply(names(meth_lab), function(nm)
  do.call(rbind, lapply(results, function(r)
    c(sigma = r$sigma, det_row(r[[nm]]))))))
rownames(det_tab) <- rep(meth_lab, each = length(results))
colnames(det_tab) <- c("sigma", paste0("x=", true_disc), "FP/rep", "FP rate")

cat("\nDetection rate within +/-", DETECT_TOL, "of each true break, and false\n")
cat("positives among the", round(n_out_mean, 1),
    "candidate midpoints outside the windows.\n")
print(as.data.frame(det_tab), digits = 3L)

det_file <- file.path(script_dir, "tab_detection.tex")

# Single-digit rates are padded with \phantom{0}, not with a space: math mode
# discards the space, and the column would lose its alignment exactly where a
# detection rate has collapsed and the reader is most likely to be comparing.
pct <- function(p) {
  s <- sprintf("%.0f", 100 * p)
  paste0("$", if (nchar(s) < 2L) "\\phantom{0}" else "", s, "\\%$")
}

det_body <- unlist(lapply(seq_along(meth_lab), function(j) {
  nm   <- names(meth_lab)[j]
  rows <- vapply(seq_along(results), function(i) {
    z <- det_row(results[[i]][[nm]])
    sprintf("  & $%.2f$ & %s & %s & %s & $%.2f$ & $%.2f\\%%$ \\\\",
            results[[i]]$sigma,
            pct(z[1]), pct(z[2]), pct(z[3]),
            z[["fp"]], 100 * z[["fp_rate"]])
  }, character(1L))
  rows[1L] <- paste0("\\multirow{", length(results), "}{*}{", meth_lab[j], "}",
                     rows[1L])
  if (j < length(meth_lab)) c(rows, "\\addlinespace") else rows
}))

writeLines(c(
  "% Generated by inst/scripts/section3_comparison.R -- do not edit by hand.",
  "% Requires \\usepackage{multirow} and \\usepackage{booktabs}.",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Discontinuity detection at each true break and false positives,",
  "by estimator and noise level. Detection rate: fraction of the",
  sprintf("$%s$ replications in which a discontinuity is flagged within",
          formatC(n_rep, big.mark = "{,}", format = "d")),
  sprintf("$\\pm%.2f$ of the true break. False positives: mean number, per",
          DETECT_TOL),
  "replication, of flagged locations falling outside all three windows, and",
  sprintf("the same figure as a percentage of the $%.1f$ candidate locations",
          n_out_mean),
  "that lie outside them on average. Discontinuities are declared at the",
  sprintf("$N-1 = %d$ midpoints between consecutive design points, which is the",
          N - 1L),
  "denominator these rates refer to.}",
  "\\label{tab:detection}",
  "\\begin{tabular}{llccccc}",
  "\\toprule",
  paste("& & \\multicolumn{3}{c}{Detection rate} &",
        "\\multicolumn{2}{c}{False positives} \\\\"),
  "\\cmidrule(lr){3-5} \\cmidrule(lr){6-7}",
  paste0("Estimator & $\\sigma$ & ",
         paste(sprintf("$x=%.1f$", true_disc), collapse = " & "),
         " & per rep. & rate \\\\"),
  "\\midrule",
  det_body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
), det_file)
cat("Saved", basename(det_file), "\n")

cat("\nMonte Carlo standard error of the per-replication RMSE:\n")
print(as.data.frame(do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$rmse_se)))), digits = 3L, row.names = FALSE)

cat("\nMedian of the per-replication RMSE (robust to the rare divergent fit):\n")
print(as.data.frame(do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$rmse_med)))), digits = 4L, row.names = FALSE)

cat("\nMean number of declared discontinuities:\n")
print(as.data.frame(do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$n_disc)))), digits = 3L, row.names = FALSE)

cat("\nCPU seconds per fit (mean over replications; user + system):\n")
print(as.data.frame(do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$time_mean)))), digits = 3L, row.names = FALSE)

cat("\nCPU seconds per fit (median, less sensitive to the odd slow fit):\n")
print(as.data.frame(do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$time_med)))), digits = 3L, row.names = FALSE)

cat("\nCPU seconds per fit, interquartile range:\n")
for (r in results)
  cat(sprintf("  sigma = %.2f : %s\n", r$sigma,
              paste(sprintf("%s [%.2f, %.2f]", colnames(r$time_iqr),
                            r$time_iqr[1L, ], r$time_iqr[2L, ]),
                    collapse = "   ")))

cat("\nFailed fits (out of", n_rep, "replications):\n")
print(as.data.frame(do.call(rbind, lapply(results, function(r)
  c(sigma = r$sigma, r$n_fail)))), row.names = FALSE)

saveRDS(results, file.path(script_dir, "section3_results.rds"))
cat("Results saved to section3_results.rds\n")

# ====================================================================
# Plotting utilities
# ====================================================================

# Fitted-curve panel: true signal + MC band + mean estimate + jump marks
draw_fit_panel <- function(res, method, col_line,
                            main = "", ylab = "", show_xlab = FALSE) {
  mn       <- res[[paste0("mean_", method)]]
  q        <- res[[paste0("q_", method)]]
  col_band <- adjustcolor(col_line, alpha.f = 0.20)
  ylim     <- range(c(mn, q, true_line), na.rm = TRUE)
  ylim     <- ylim + c(-1, 1) * 0.06 * diff(ylim)

  plot(NA, xlim = c(0, 1), ylim = ylim,
       xlab = if (show_xlab) "x" else "",
       ylab = ylab, main = main, xaxs = "i", las = 1L)
  lines(x_grid, true_line, lty = 2L, col = "grey40", lwd = 0.9)
  polygon(c(x_grid, rev(x_grid)), c(q[1L, ], rev(q[2L, ])),
          col = col_band, border = NA)
  lines(x_grid, mn, col = col_line, lwd = 1.8)
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
  plot(NA, xlim = c(0, 1), ylim = c(0, n_rep),
       xaxs = "i", yaxs = "i",
       xlab = if (show_xlab) "x" else "", ylab = "",
       main = "", axes = FALSE)
  rect(brk[-length(brk)], 0L, brk[-1L], cnt,
       col = adjustcolor(col_line, alpha.f = 0.65), border = NA)
  axis(1L, at = c(0, 0.3, 0.4, 0.6, 1))
  axis(2L, las = 1L, at = pretty(c(0L, n_rep), n = 3L))
  abline(v = true_disc, lty = 3L, col = "grey20")
}

# ====================================================================
# Figures: one 2 x 3 PDF per method
# ====================================================================

methods   <- c("cssd",      "atw",                        "mle",         "ss")
meth_labs <- list(
  "CSSD (CV)",
  quote("ATW (CV-" * lambda * ")"),
  "MLE",
  "Smooth. spline"
)
col_lines <- c("steelblue", "darkgreen", "darkorange3", "purple4")

for (m in seq_along(methods)) {
  meth     <- methods[m]
  col_line <- col_lines[m]
  lab      <- meth_labs[[m]]
  fname    <- paste0("fig_", meth, ".pdf")
  has_disc <- meth != "ss"

  pdf(file.path(script_dir, fname), width = 7,
      height = if (has_disc) 4.5 else 2.5)
  if (has_disc) {
    layout(matrix(seq_len(6L), nrow = 2L, ncol = 3L, byrow = TRUE),
           heights = c(3.5, 1.8))
    op <- par(mar = c(1.5, 3.5, 2.5, 0.5), mgp = c(2, 0.5, 0),
              tcl = -0.25, cex = 0.85, cex.main = 0.95)
  } else {
    op <- par(mfrow = c(1L, 3L), mar = c(3.5, 3.5, 2.5, 0.5),
              mgp = c(2, 0.5, 0), tcl = -0.25, cex = 0.85, cex.main = 0.95)
  }

  for (s in seq_along(sigmas)) {
    draw_fit_panel(results[[s]], method = meth, col_line = col_line,
                   main      = bquote(.(lab) * "  (" * sigma * " = " * .(sigmas[s]) * ")"),
                   ylab      = if (s == 1L) "y" else "",
                   show_xlab = !has_disc)
  }

  if (has_disc) {
    par(mar = c(3, 3.5, 0.3, 0.5))
    for (s in seq_along(sigmas)) {
      draw_hist_panel(results[[s]][[paste0("disc_", meth)]],
                      col_line  = col_line,
                      show_xlab = TRUE)
      if (s == 1L) mtext("count", side = 2L, line = 2.5, cex = 0.72)
    }
  }

  par(op)
  dev.off()
  cat("Saved", fname, "\n")
}

# ====================================================================
# Combined figure: the three discontinuity-aware methods only.
# 6 rows x 3 cols: one curve/histogram pair per method.
#
# The classical smoothing spline is deliberately absent.  It cannot declare a
# discontinuity, so it would contribute an odd row out with no histogram
# beneath it, and its role here is as a benchmark for the RMSE table rather
# than as a competitor in locating the jumps.  It keeps its own figure
# (fig_ss.pdf) and its column in every table.
# ====================================================================

fig_methods <- setdiff(methods, "ss")

pdf(file.path(script_dir, "fig_combined.pdf"), width = 7, height = 7.6)
mat <- matrix(seq_len(6L * 3L), nrow = 6L, ncol = 3L, byrow = TRUE)
layout(mat, heights = rep(c(3.5, 1.8), length(fig_methods)))
op <- par(mgp = c(2, 0.5, 0), tcl = -0.25, cex = 0.80, cex.main = 0.93)

for (meth in fig_methods) {
  m        <- match(meth, methods)
  col_line <- col_lines[m]
  lab      <- meth_labs[[m]]

  par(mar = c(1.5, 3.5, 2.5, 0.5))
  for (s in seq_along(sigmas)) {
    draw_fit_panel(results[[s]], method = meth, col_line = col_line,
                   main      = bquote(.(lab) * "  (" * sigma * " = " * .(sigmas[s]) * ")"),
                   ylab      = if (s == 1L) "y" else "",
                   show_xlab = FALSE)
  }

  par(mar = c(3, 3.5, 0.3, 0.5))
  for (s in seq_along(sigmas)) {
    draw_hist_panel(results[[s]][[paste0("disc_", meth)]],
                    col_line  = col_line,
                    show_xlab = TRUE)
    if (s == 1L) mtext("count", side = 2L, line = 2.5, cex = 0.70)
  }
}

par(op)
invisible(dev.off())
cat("Saved fig_combined.pdf\n")

# ====================================================================
# Figure: pointwise Bias and MSE (2 rows x 3 columns, all sigma levels)
# ====================================================================

sigma_ltys <- c(1L, 2L, 4L)

bias_list <- lapply(seq_along(methods), function(m)
  lapply(results, function(res)
    res[[paste0("mean_", methods[m])]] - true_line))

mse_list <- lapply(seq_along(methods), function(m)
  lapply(results, function(res) res[[paste0("mse_", methods[m])]]))

bias_ylim <- local({
  r <- range(unlist(bias_list), na.rm = TRUE)
  max(abs(r)) * c(-1.1, 1.1)
})
mse_ylim  <- c(0, max(unlist(mse_list), na.rm = TRUE) * 1.1)

sig_labs  <- as.expression(lapply(sigmas, function(s) bquote(sigma == .(s))))

pdf(file.path(script_dir, "fig_bias_mse.pdf"), width = 8, height = 4)
layout(matrix(seq_len(8L), nrow = 2L, ncol = 4L, byrow = TRUE))
op <- par(mgp = c(2.4, 0.55, 0), tcl = -0.25, cex = 0.85, cex.main = 0.95)

# Row 1: Bias
par(mar = c(1.5, 3.8, 2.5, 0.5))
for (m in seq_along(methods)) {
  lab      <- meth_labs[[m]]
  plot(NA, xlim = c(0, 1), ylim = bias_ylim,
       xlab = "", ylab = if (m == 1L) "Bias" else "",
       main = bquote(.(lab)), xaxs = "i", las = 1L)
  abline(h = 0, col = "grey70", lwd = 0.8)
  abline(v = true_disc, lty = 3L, col = "grey20")
  for (s in seq_along(sigmas)) {
    col_line <- col_lines[s]
    lines(x_grid, bias_list[[m]][[s]], col = col_line,
            lty = 1, lwd = 1.5)
  }
  if (m == 1L)
    legend("bottomright", legend = sig_labs, lty = 1,
           col = col_lines, lwd = 1.5, bty = "n", cex = 0.78)
}

# Row 2: MSE
par(mar = c(3.5, 3.8, 0.8, 0.5))
for (m in seq_along(methods)) {
  plot(NA, xlim = c(0, 1), ylim = mse_ylim,
       xlab = "x", ylab = if (m == 1L) "MSE" else "",
       xaxs = "i", las = 1L, main = "")
  abline(v = true_disc, lty = 3L, col = "grey20")
  for (s in seq_along(sigmas)) {
    col_line <- col_lines[s]
    lines(x_grid, mse_list[[m]][[s]], col = col_line,
          lty = 1, lwd = 1.5)
  }
}

par(op)
invisible(dev.off())
cat("Saved fig_bias_mse.pdf\n")

# ====================================================================
# Figure: computing time per fit
#
# One box per method and noise level, on a log scale because the three differ
# by more than an order of magnitude and the smoothing spline is far below the
# rest.  CPU seconds, so the boxes do not depend on how many replications were
# running concurrently.  Drawn from time_all, which keeps every measurement.
# ====================================================================

time_meths <- c("cssd", "mle", "atw")
time_labs  <- c("CSSD (CV)", "MLE", expression("ATW (CV-" * lambda * ")"))
time_cols  <- col_lines[match(time_meths, methods)]

pdf(file.path(script_dir, "fig_timing.pdf"), width = 7, height = 3.2)
op <- par(mfrow = c(1L, length(sigmas)), mar = c(3.2, 3.8, 2.2, 0.6),
          mgp = c(2.4, 0.6, 0), tcl = -0.25, cex = 0.85, cex.main = 0.95)

ylim <- range(unlist(lapply(results, function(r) r$time_all[, time_meths])),
              na.rm = TRUE)
ylim[1L] <- max(ylim[1L], 1e-3)          # keep the log scale finite

for (s in seq_along(sigmas)) {
  dat <- lapply(time_meths, function(m) {
    v <- results[[s]]$time_all[, m]
    v[is.finite(v) & v > 0]
  })
  boxplot(dat, log = "y", ylim = ylim, names = rep("", length(time_meths)),
          col = adjustcolor(time_cols, alpha.f = 0.35), border = time_cols,
          ylab = if (s == 1L) "CPU seconds per fit" else "",
          main = bquote(sigma == .(sigmas[s])), las = 1L, outcex = 0.4)
  axis(1L, at = seq_along(time_meths), labels = FALSE)
  mtext(time_labs, side = 1L, at = seq_along(time_meths), line = 0.9, cex = 0.72)
}

par(op)
invisible(dev.off())
cat("Saved fig_timing.pdf
")
