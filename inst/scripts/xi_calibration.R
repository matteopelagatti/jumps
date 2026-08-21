## ===========================================================================
##  Calibration of the extended-BIC parameter xi, and comparison against
##  cross-validation, for the domain-warping smoothing splines with jumps.
##
##  Companion to Section 3.  Produces:
##    * risk curves in xi (RMSE against the true signal, and balanced accuracy
##      of jump detection) with PAIRED Monte Carlo standard errors;
##    * the one-standard-error interval for the optimal xi in every design
##      cell, and the minimax-regret xi across cells;
##    * for ssj_atw, a comparison of xi-calibrated EBIC against 5-fold CV over
##      the budget M, over (lambda, M) jointly, and a hybrid that takes lambda
##      from CV and M from EBIC;
##    * the regression of the effective degrees of freedom on the number of
##      active warps, which quantifies how much of the search cost the plain
##      BIC term already charges (needed to reconcile the calibrated xi with
##      the Chen and Chen, 2008, consistency floor).
##
##  DESIGN NOTES
##
##  1. xi enters ONLY the selection among grid points, never a fit.  The
##     regularisation path is therefore computed once per replication and every
##     xi is evaluated by re-selecting from the stored path.  The whole
##     xi-curve costs one simulation, and all xi are compared on IDENTICAL fits
##     -- a fully paired design, so the standard error of a DIFFERENCE between
##     two xi is far smaller than the standard error of either level.
##
##  2. The calibration signal is NOT the Bessel signal used for evaluation in
##     Section 3: jump locations and sizes are redrawn at every replication.
##     The chosen xi is then validated out-of-sample on the Section 3 design.
##
##  3. Common random numbers: the standard normal draws are generated before
##     the noise level is applied, so every sigma sees the same design and the
##     same standardised errors, and the comparison across sigma is paired too.
##
##  4. Two different jump counts appear below, on purpose.
##       k_pkg  is the count the package itself puts in the ebic penalty (the
##              raw number of warps above its quasi-zero).  It is recovered
##              EXACTLY from the returned criteria, as
##                 k_pkg = (ebic - bic) / (2 * xi_default * log(n-1)),
##              so the xi calibrated here transfers verbatim to auto_ssj_*().
##       k_decl is the number of jumps DECLARED by the package's reporting rule
##              ssj_breaks() (a warp must clear 5% of the largest warp and 30%
##              of the mean spacing).  This is the paper's rule and is the one
##              that belongs in the detection metric.
##     Calibrating the penalty on k_pkg while scoring detection on k_decl is
##     what makes the answer usable without re-engineering the estimators.
##
##  5. ssj_fix is deliberately absent.  It has no ebic among its information
##     criteria and parameterises the warps as additional standard deviations
##     rather than as a budget on gamma, so there is no xi to calibrate.
##
##  6. Cross-validation arms are computed for ssj_atw only.  Its f-step is an
##     O(n) banded solve, so K folds over a two-dimensional grid is affordable;
##     the same design over ssj_mle or ssj_kfs would mean K * |lambda| * |M|
##     constrained optimisations per replication, which is not.  lambda is
##     shared across estimators, as in the paper, so the CV-selected lambda
##     reaches all of them at the cost of the ATW folds alone.
##
##  USAGE
##    From a shell:
##       Rscript xi_calibration.R probe     # time it and extrapolate
##       Rscript xi_calibration.R run       # checkpoints after every cell
##       Rscript xi_calibration.R analyse   # tables from the checkpoints
##    From RStudio, set MODE first and source the file:
##       MODE <- "run"; source("inst/scripts/xi_calibration.R")
##    (an existing MODE always wins over the command line, so re-sourcing in
##    the same session keeps doing what you asked for).
## ===========================================================================

suppressMessages(library(jumps))
library(parallel)

## ---------------------------------------------------------------- CONFIG ---

if (!exists("MODE")) {
  a <- commandArgs(trailingOnly = TRUE)
  MODE <- if (length(a) && nzchar(a[1])) a[1] else "probe"
}

OUT_DIR <- file.path(getwd(), "xi_calibration_out")

## Replications are dispatched in this many chunks per cell, purely so that
## progress appears while a cell is running; the cost of the extra round trips
## is negligible against the fits themselves.
PROGRESS_CHUNKS <- 8L

ESTIMATORS <- c("mle", "kfs", "atw")   # the three that have an ebic
N_GRID     <- c(50L, 100L, 200L, 400L)
SIGMA_GRID <- c(0.05, 0.10, 0.15)
J_GRID     <- c(0L, 1L, 3L, 5L)        # J = 0 measures false discoveries

## replications per cell; large n costs much more, so allow a taper
R_BY_N <- c("50" = 500L, "100" = 500L, "200" = 300L, "400" = 150L)

XI_GRID <- seq(0, 1, by = 0.05)
XI_REF  <- 0.25       # used ONLY to pick the lambda that the fixed-lambda
                      # estimators inherit from ssj_mle, and as the EBIC arm in
                      # the CV comparison.  Held fixed so the calibration of xi
                      # is not circular.

## The package's own ebic_xi defaults, needed to recover k_pkg exactly.  Read
## from the functions rather than hard-coded: the defaults are themselves an
## output of this experiment, so a literal here would silently go stale the
## moment they were updated, and k_pkg would come back scaled by the ratio.
XI_DEFAULT <- vapply(list(mle = ssj_mle, kfs = ssj_kfs, atw = ssj_atw),
                     function(f) eval(formals(f)$ebic_xi), numeric(1))
stopifnot(all(is.finite(XI_DEFAULT)), all(XI_DEFAULT > 0))

## Budget grid, in units of the range of x (equal to 1 here by construction).
## It runs to twice the range: the pilot pinned the selection to the top of a
## grid that stopped at the range itself for small xi, and a censored selection
## is not an answer.  Boundary hits are recorded and reported.
M_GRID <- seq(0, 2, length.out = 25L)

LAMBDA_SOURCES <- c("ml", "cv")        # for kfs and atw; mle estimates its own
LAMBDA_GRID    <- 10^seq(-8, -1, length.out = 8L)
K_FOLD         <- 5L

N_CORES <- max(1L, detectCores() - 1L)

## Declaration rule for a jump: the package's own, i.e. the paper's rule.
BREAK_ARGS <- list(frac_total = 0.05, frac_delta = 0.30)

X_EVAL <- seq(0, 1, length.out = 500L)

## --- analysis-only knobs -------------------------------------------------
## Trim used when averaging risk across replications.  0 reproduces the
## plain mean; 0.02 makes every table robust to the rare catastrophic fit
## that the extreme-fit diagnostic reports (see section 2 of MODE=analyse).
ROBUST_TRIM <- 0
## Print the 48-line per-cell table for each arm, or only the summaries.
SHOW_CELLS  <- TRUE

## ------------------------------------------------------------- DATA GEN ---

## Smooth background, deliberately unlike besselJ(20x, 1).
background <- function(x) sin(2*pi*x) + 0.5*sin(5*pi*x)

make_data <- function(n, J, rep, sigma) {
  set.seed(20260816L + 1000L*J + rep)
  x <- sort(runif(n))
  z <- rnorm(n)                       # drawn before sigma is applied: CRN
  if (J > 0L) {
    locs  <- sort(runif(J, 0.12, 0.88))
    sizes <- sample(c(-1, 1), J, TRUE) * runif(J, 0.30, 0.70)
    f <- function(u) background(u) + as.numeric(outer(u, locs, ">=") %*% sizes)
  } else {
    locs <- numeric(0); f <- background
  }
  list(x = x, y = f(x) + sigma*z, f = f, locs = locs)
}

## ------------------------------------------------------------- FITTING ----

## The regularisation path over M_GRID, warm started from the previous point.
## Entries are NULL where the fit failed.
fit_path <- function(est, y, x, lambda, M_grid = M_GRID) {
  out  <- vector("list", length(M_grid))
  init <- NULL
  for (j in seq_along(M_grid)) {
    f <- tryCatch(switch(est,
        mle = ssj_mle(y, x, maxsum = M_grid[j], parinit = init),
        kfs = ssj_kfs(y, x, lambda = lambda, maxsum = M_grid[j], parinit = init),
        atw = ssj_atw(y, x, lambda = lambda, maxsum = M_grid[j],
                      gamma_init = if (is.null(init)) numeric(0) else init)),
      error = function(e) NULL)
    if (is.null(f)) { init <- NULL; next }
    ## warm start, in each estimator's own parameter layout:
    ##   ssj_mle : (sigma, lambda, gamma_1..gamma_{n-1})   length n+1
    ##   ssj_kfs : (sigma, gamma_1..gamma_{n-1})           length n
    ##   ssj_atw : gamma only
    init <- switch(est,
      mle = c(unname(f$pars["sigma"]), unname(f$pars["lambda"]), f$gamma),
      kfs = c(unname(f$pars["sigma"]), f$gamma),
      atw = f$gamma)
    out[[j]] <- f
  }
  out
}

## The jump count that the package charged in its own ebic, recovered exactly.
k_package <- function(f, est) {
  nn <- f$nobs
  d  <- unname(f$ic["ebic"] - f$ic["bic"])
  if (!is.finite(d) || nn < 3L) return(NA_real_)
  d / (2 * XI_DEFAULT[[est]] * log(nn - 1))
}

## ------------------------------------------------------------- SCORING ----

## ssj_breaks() dates a break at the index of the point AFTER it, so a jump
## inside (x_i, x_{i+1}) is break index i+1.
true_break_index <- function(locs, xf) {
  if (!length(locs)) return(integer(0))
  vapply(locs, function(L)
    max(1L, min(length(xf) - 1L, findInterval(L, xf))), integer(1)) + 1L
}

score_fit <- function(f, truth, true_idx, ncand, M_j) {
  xe <- pmin(pmax(X_EVAL, min(f$x)), max(f$x))
  p  <- tryCatch(predict(f, xe), error = function(e) rep(NA_real_, length(xe)))
  br <- tryCatch(do.call(ssj_breaks, c(list(f), BREAK_ARGS))$index,
                 error = function(e) integer(0))
  pos_window <- unique(unlist(lapply(true_idx, function(t) (t-1L):(t+1L))))
  neg <- setdiff(seq.int(2L, ncand + 1L), pos_window)
  rec  <- if (length(true_idx))
            mean(vapply(true_idx, function(t) any(abs(br - t) <= 1L), TRUE))
          else NA_real_
  spec <- if (length(neg)) mean(!(neg %in% br)) else NA_real_
  c(rmse = sqrt(mean((p - truth)^2)),
    ba   = if (is.na(rec)) spec else (rec + spec)/2,
    rec  = rec, spec = spec, k_decl = length(br), M = M_j)
}

## ------------------------------------------------------- CV MACHINERY -----

## K-fold CV surface for ATW over (lambda, M), same folds throughout so the
## comparison is paired.  Returns a |lambda| x |M| matrix of mean squared
## prediction error.
cv_surface_atw <- function(y, x, fold, lam_grid = LAMBDA_GRID, M_grid = M_GRID) {
  err <- cnt <- matrix(0, length(lam_grid), length(M_grid))
  for (k in sort(unique(fold))) {
    tr <- fold != k; te <- !tr
    if (!any(te) || sum(tr) < 5L) next
    for (a in seq_along(lam_grid)) {
      pk <- fit_path("atw", y[tr], x[tr], lam_grid[a], M_grid)
      for (j in seq_along(M_grid)) {
        f <- pk[[j]]; if (is.null(f)) next
        xe <- pmin(pmax(x[te], min(f$x)), max(f$x))
        pr <- tryCatch(predict(f, xe), error = function(e) NULL)
        if (is.null(pr) || anyNA(pr)) next
        err[a, j] <- err[a, j] + sum((y[te] - pr)^2)
        cnt[a, j] <- cnt[a, j] + sum(te)
      }
    }
  }
  ifelse(cnt > 0, err/cnt, Inf)
}

## ------------------------------------------------------- ONE REPLICATION --

run_rep <- function(rep, n, J, sigma) {
  d <- make_data(n, J, rep, sigma)
  y <- d$y; x <- d$x
  truth <- d$f(X_EVAL)
  rows <- list(); dfk <- list()

  add <- function(est, lsrc, sel, xi, s, kp) {
    rows[[length(rows) + 1L]] <<- data.frame(
      n = n, sigma = sigma, J = J, rep = rep, estimator = est,
      lambda_src = lsrc, selector = sel, xi = xi,
      rmse = unname(s["rmse"]), ba = unname(s["ba"]), rec = unname(s["rec"]),
      spec = unname(s["spec"]), k_decl = unname(s["k_decl"]), k_pkg = kp,
      M = unname(s["M"]),
      at_boundary = unname(s["M"]) >= M_GRID[length(M_GRID)] - 1e-12,
      stringsAsFactors = FALSE)
  }

  ## Evaluate every xi on a stored path.  Returns the chosen index per xi.
  pick <- function(bic, kp, nn) {
    vapply(XI_GRID, function(xi) which.min(bic + 2*xi*kp*log(nn - 1)), integer(1))
  }

  harvest <- function(est, lsrc, pp) {
    ok <- !vapply(pp, is.null, TRUE)
    if (sum(ok) < 3L) return(NULL)
    ff <- pp[ok]; Ms <- M_GRID[ok]
    xf <- ff[[1]]$x
    tix <- true_break_index(d$locs, xf); ncand <- length(xf) - 1L
    nn  <- ff[[1]]$nobs
    bic <- vapply(ff, function(f) unname(f$ic["bic"]), numeric(1))
    kp  <- vapply(ff, k_package, numeric(1), est = est)
    dv  <- vapply(ff, function(f) f$df, numeric(1))
    if (anyNA(kp) || anyNA(bic)) return(NULL)
    sc  <- lapply(seq_along(ff), function(i) score_fit(ff[[i]], truth, tix, ncand, Ms[i]))
    dfk[[length(dfk) + 1L]] <<- data.frame(
      n = n, sigma = sigma, J = J, rep = rep, estimator = est,
      lambda_src = lsrc, df = dv, k_pkg = kp,
      k_decl = vapply(sc, `[[`, numeric(1), "k_decl"), stringsAsFactors = FALSE)
    idx <- pick(bic, kp, nn)
    for (i in seq_along(XI_GRID)) add(est, lsrc, "ebic", XI_GRID[i], sc[[idx[i]]], kp[idx[i]])
    ib <- which.min(vapply(sc, `[[`, numeric(1), "rmse"))
    add(est, lsrc, "oracle", NA_real_, sc[[ib]], kp[ib])
    list(ok = ok, ff = ff, sc = sc, bic = bic, kp = kp, nn = nn)
  }

  ## ---- ssj_mle: estimates lambda itself -----------------------------------
  h_mle <- harvest("mle", "self", fit_path("mle", y, x, NA_real_))
  if (is.null(h_mle)) return(NULL)

  ## the lambda that the fixed-lambda estimators inherit, chosen at XI_REF
  b_ref  <- which.min(h_mle$bic + 2*XI_REF*h_mle$kp*log(h_mle$nn - 1))
  lam_ml <- unname(h_mle$ff[[b_ref]]$pars["lambda"])

  ## ---- CV surface for ATW: yields lambda_cv and the CV arms ---------------
  fold  <- sample(rep_len(seq_len(K_FOLD), n))
  cvs   <- tryCatch(cv_surface_atw(y, x, fold), error = function(e) NULL)
  if (is.null(cvs) || all(!is.finite(cvs))) return(NULL)
  b2     <- which(cvs == min(cvs), arr.ind = TRUE)[1, ]
  lam_cv <- LAMBDA_GRID[b2[["row"]]]
  j_cv2  <- b2[["col"]]                                  # M chosen jointly
  a_ml   <- which.min(abs(log(LAMBDA_GRID) - log(lam_ml)))
  j_cv1  <- which.min(cvs[a_ml, ])                       # M at the ML lambda

  lams <- c(ml = lam_ml, cv = lam_cv)[LAMBDA_SOURCES]

  ## ---- ssj_kfs and ssj_atw at each lambda source --------------------------
  for (lsrc in names(lams)) {
    lam <- lams[[lsrc]]
    if (!is.finite(lam) || lam <= 0) next
    for (est in intersect(ESTIMATORS, c("kfs", "atw"))) {
      h <- harvest(est, lsrc, fit_path(est, y, x, lam))
      if (is.null(h) || est != "atw") next

      ## CV arms (ATW only).  Map the CV-chosen grid position onto the
      ## positions that actually fitted.
      at <- function(j) { p <- which.min(abs(which(h$ok) - j)); h$sc[[p]] }
      kp_at <- function(j) { p <- which.min(abs(which(h$ok) - j)); h$kp[p] }
      if (lsrc == "ml") {
        add(est, lsrc, "cv_M", NA_real_, at(j_cv1), kp_at(j_cv1))
      } else {
        add(est, lsrc, "cv_lambda_M", NA_real_, at(j_cv2), kp_at(j_cv2))
        b <- which.min(h$bic + 2*XI_REF*h$kp*log(h$nn - 1))
        add(est, lsrc, "cv_lambda_ebic_M", XI_REF, h$sc[[b]], h$kp[b])
      }
    }
  }

  list(rows = do.call(rbind, rows), dfk = do.call(rbind, dfk))
}

## ------------------------------------------------------------- DRIVER -----

cells <- expand.grid(n = N_GRID, sigma = SIGMA_GRID, J = J_GRID,
                     KEEP.OUT.ATTRS = FALSE)
cells <- cells[order(cells$n, cells$J, cells$sigma), ]

EXPORTS <- c("ESTIMATORS","XI_GRID","XI_REF","XI_DEFAULT","M_GRID","LAMBDA_GRID",
             "LAMBDA_SOURCES","K_FOLD","BREAK_ARGS","X_EVAL","background",
             "make_data","fit_path","k_package","true_break_index","score_fit",
             "cv_surface_atw","run_rep")

if (MODE == "probe") {

  cat("Timing probe (serial, 2 replications each).\n\n")
  probe <- vapply(c(min(N_GRID), max(N_GRID)), function(nn) {
    t0 <- Sys.time()
    invisible(lapply(1:2, function(r)
      run_rep(r, nn, max(J_GRID), stats::median(SIGMA_GRID))))
    per <- as.numeric(difftime(Sys.time(), t0, units = "secs"))/2
    cat(sprintf("  n = %3d, J = %d : %7.2f s per replication\n",
                nn, max(J_GRID), per))
    per
  }, numeric(1))
  n_lo <- min(N_GRID); n_hi <- max(N_GRID)
  cost <- function(nn) if (n_hi == n_lo) probe[1] else
    probe[1] * (probe[2]/probe[1])^((log(nn) - log(n_lo))/(log(n_hi) - log(n_lo)))
  tot <- sum(vapply(seq_len(nrow(cells)), function(i)
    cost(cells$n[i]) * R_BY_N[[as.character(cells$n[i])]], numeric(1)))
  cat(sprintf("\n  %d cells, %s replications in total\n", nrow(cells),
              format(sum(R_BY_N[as.character(cells$n)]), big.mark = ",")))
  cat(sprintf("  estimated serial   : %6.1f h\n", tot/3600))
  cat(sprintf("  estimated %2d cores : %6.1f h  (perfect scaling; expect 1.5-2x more)\n",
              N_CORES, tot/3600/N_CORES))
  cat("\n  Tune R_BY_N / N_GRID / J_GRID, then rerun with argument \"run\".\n")

} else if (MODE == "run") {

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  cl <- makeCluster(N_CORES)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  invisible(clusterEvalQ(cl, suppressMessages(library(jumps))))
  clusterExport(cl, EXPORTS, envir = globalenv())

  hms <- function(s) sprintf("%d:%02d", as.integer(s) %/% 3600,
                             (as.integer(s) %% 3600) %/% 60)
  say <- function(...) { cat(...); flush.console() }

  t_start   <- Sys.time()
  work_done <- 0                       # replications finished this session
  work_todo <- sum(vapply(seq_len(nrow(cells)), function(i) {
    tag <- sprintf("n%03d_s%03.0f_J%d", cells$n[i], 100*cells$sigma[i], cells$J[i])
    if (file.exists(file.path(OUT_DIR, paste0("cell_", tag, ".rds")))) 0
    else R_BY_N[[as.character(cells$n[i])]] }, numeric(1)))
  say(sprintf("%d cells, %s replications still to do.\n\n",
              nrow(cells), format(work_todo, big.mark = ",")))

  for (i in seq_len(nrow(cells))) {
    nn <- cells$n[i]; sg <- cells$sigma[i]; jj <- cells$J[i]
    tag <- sprintf("n%03d_s%03.0f_J%d", nn, 100*sg, jj)
    fn  <- file.path(OUT_DIR, paste0("cell_", tag, ".rds"))
    if (file.exists(fn)) { say(sprintf("[%2d/%2d] %s  skip (done)\n",
                                       i, nrow(cells), tag)); next }
    R  <- R_BY_N[[as.character(nn)]]
    t0 <- Sys.time()
    say(sprintf("[%2d/%2d] %s  n=%3d sigma=%.2f J=%d  R=%d  ",
                i, nrow(cells), tag, nn, sg, jj, R))

    ## dispatched in chunks so that a dot appears every few per cent
    idx <- split(seq_len(R), cut(seq_len(R), PROGRESS_CHUNKS, labels = FALSE))
    rr  <- list()
    for (ch in idx) {
      rr <- c(rr, parLapplyLB(cl, ch, run_rep, n = nn, J = jj, sigma = sg))
      say(".")
    }

    bad  <- vapply(rr, is.null, TRUE)
    secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    saveRDS(list(rows  = do.call(rbind, lapply(rr[!bad], `[[`, "rows")),
                 dfk   = do.call(rbind, lapply(rr[!bad], `[[`, "dfk")),
                 nfail = sum(bad), R = R, secs = secs), fn)

    work_done <- work_done + R
    elapsed   <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    eta       <- elapsed/work_done * (work_todo - work_done)
    say(sprintf("  %d/%d ok, %4.0fs | elapsed %s | ETA %s\n",
                R - sum(bad), R, secs, hms(elapsed), hms(eta)))
  }
  stopCluster(cl); on.exit()
  cat(sprintf("\nDone in %s. Now run with MODE <- \"analyse\".\n",
              hms(as.numeric(difftime(Sys.time(), t_start, units = "secs")))))

} else if (MODE == "analyse") {

  fs <- list.files(OUT_DIR, "^cell_.*\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no checkpoints in ", OUT_DIR, "; run with \"run\" first")
  parts <- lapply(fs, readRDS)
  D  <- do.call(rbind, lapply(parts, `[[`, "rows"))
  FK <- do.call(rbind, lapply(parts, `[[`, "dfk"))
  cat(sprintf("%s rows from %d cells; %d replications failed\n\n",
              format(nrow(D), big.mark = ","), length(fs),
              sum(vapply(parts, `[[`, numeric(1), "nfail"))))

  cellid <- function(z) paste(z$n, z$sigma, z$J, sep = "/")

  ## ======================================================================
  ##  1. Detection metrics that survive the class imbalance
  ##
  ##  BA = (recall + specificity)/2 does not work here.  The candidate space
  ##  has n-1 members of which only about 3J lie near a true jump, so the
  ##  negative class outnumbers the positive one by two orders of magnitude:
  ##  a handful of false positives moves specificity by O(1/n), and BA
  ##  collapses to recall/2 + 1/2.  It is therefore maximised by the smallest
  ##  penalty when J >= 1, and by the largest when J = 0 (where it degenerates
  ##  to pure specificity) -- two artefacts of the metric, not two findings.
  ##
  ##  The Matthews correlation coefficient uses all four cells of the
  ##  confusion matrix and is the standard choice under this imbalance.
  ##  Everything it needs is reconstructed from what was already stored, so no
  ##  refit is required:  TP = recall * J,  FN = J - TP,  and with
  ##  nneg = (n-1) - 3J negatives,  TN = specificity * nneg  and
  ##  FP = (1 - specificity) * nneg.  MCC is undefined at J = 0 (there are no
  ##  positives), which is exactly right: those cells are reported separately
  ##  as false-discovery counts instead of being allowed to vote.
  ## ======================================================================
  add_detection <- function(z) {
    nneg <- pmax((z$n - 1) - 3*z$J, 1)
    TP <- ifelse(z$J > 0, z$rec * z$J, 0)
    FN <- pmax(z$J - TP, 0)
    FP <- (1 - z$spec) * nneg
    TN <- z$spec * nneg
    den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    z$mcc   <- ifelse(z$J == 0 | !is.finite(den) | den <= 0, NA_real_,
                      (TP*TN - FP*FN) / den)
    z$f1    <- ifelse(z$J == 0 | (2*TP + FP + FN) <= 0, NA_real_,
                      2*TP / (2*TP + FP + FN))
    z$fdisc <- FP
    z
  }
  D <- add_detection(D)

  ## ======================================================================
  ##  2. Extreme fits
  ##
  ##  EBIC selects on the likelihood, so it can pick a fit that is excellent
  ##  at the design points and wild between them: a large warp makes a long
  ##  latent interval, and the cubic overshoots when projected back to the
  ##  physical domain.  The oracle and the CV arms cannot do this, because
  ##  both select on an error criterion -- which is why only the ebic rows
  ##  show it.  One such fit moves a 150-replication mean by an order of
  ##  magnitude, so the mean is not a safe summary until these are counted.
  ## ======================================================================
  key <- with(D, paste(estimator, lambda_src, selector, n, sigma, J, sep = "|"))
  med <- ave(D$rmse, key, FUN = function(v) median(v, na.rm = TRUE))
  D$extreme <- is.finite(D$rmse) & D$rmse > 10*med

  cat("=== extreme fits (RMSE more than 10x the median of their own cell) ===\n")
  ex <- D[D$extreme, ]
  if (!nrow(ex)) {
    cat("  none\n\n")
  } else {
    g <- split(ex, paste(ex$estimator, ex$lambda_src, ex$selector,
                         ex$n, ex$sigma, ex$J, sep = " "))
    tb <- do.call(rbind, lapply(names(g), function(k) {
      z <- g[[k]]
      data.frame(cell = k, n_bad = nrow(z),
                 max_rmse = max(z$rmse, na.rm = TRUE),
                 reps = paste(utils::head(sort(unique(z$rep)), 6), collapse = ","),
                 xi = paste(utils::head(sort(unique(z$xi)), 6), collapse = ","),
                 stringsAsFactors = FALSE) }))
    print(tb[order(-tb$max_rmse), ], row.names = FALSE, digits = 4)
    cat(sprintf("\n  %d extreme rows out of %s (%.4f%%), in %d cell/selector groups.\n",
                nrow(ex), format(nrow(D), big.mark = ","),
                100*nrow(ex)/nrow(D), nrow(tb)))
    cat("  The replication ids above reproduce them: make_data(n, J, rep, sigma).\n\n")
    cat("  Affected groups, mean against robust summaries:\n")
    for (k in names(g)) {
      z <- D[paste(D$estimator, D$lambda_src, D$selector,
                   D$n, D$sigma, D$J, sep = " ") == k, ]
      cat(sprintf("    %-42s mean %9.4f | median %7.4f | 10%% trimmed %7.4f\n",
                  k, mean(z$rmse, na.rm = TRUE), median(z$rmse, na.rm = TRUE),
                  mean(z$rmse, trim = 0.10, na.rm = TRUE)))
    }
    if (ROBUST_TRIM <= 0)
      cat("\n  Set ROBUST_TRIM (top of the file) to 0.02 and re-run to make every\n",
          " table below use a trimmed mean instead.\n", sep = "")
    cat("\n")
  }

  ## ======================================================================
  ##  3. Per-cell optima, with PAIRED one-standard-error intervals
  ## ======================================================================
  one_se <- function(sub, metric, maximise = FALSE) {
    sub <- sub[, c("rep", "xi", metric)]
    if (!any(is.finite(sub[[metric]]))) return(NULL)
    w  <- reshape(sub, idvar = "rep", timevar = "xi", direction = "wide")
    m  <- as.matrix(w[, -1, drop = FALSE])
    xi <- as.numeric(sub(paste0("^", metric, "\\."), "", colnames(m)))
    o  <- order(xi); xi <- xi[o]; m <- m[, o, drop = FALSE]
    mu <- apply(m, 2, mean, trim = ROBUST_TRIM, na.rm = TRUE)
    if (all(is.na(mu))) return(NULL)
    b   <- if (maximise) which.max(mu) else which.min(mu)
    dif <- if (maximise) m[, b] - m else m - m[, b]      # paired regret
    se  <- apply(dif, 2, function(v) {
             v <- v[is.finite(v)]
             if (length(v) < 2L) NA_real_ else sd(v)/sqrt(length(v)) })
    gap <- if (maximise) mu[b] - mu else mu - mu[b]
    names(gap) <- xi
    inside <- which(gap <= se | !is.finite(se))
    list(best = xi[b], lo = min(xi[inside]), hi = max(xi[inside]), regret = gap)
  }

  minimax <- function(RG, lab) {
    if (!length(RG)) { cat(sprintf("  minimax-regret xi on %-4s : no cells\n", lab))
                       return(NA_real_) }
    M <- do.call(rbind, RG)
    ## a regret column can be entirely missing (MCC is undefined whenever no
    ## jump is declared and none is near one); such a xi cannot be compared
    ## and is dropped rather than allowed to win by default
    worst <- apply(M, 2, function(v)
      if (!any(is.finite(v))) NA_real_ else max(v, na.rm = TRUE))
    if (!any(is.finite(worst))) {
      cat(sprintf("  minimax-regret xi on %-4s : not identified (no finite regrets)\n", lab))
      return(NA_real_)
    }
    xi <- as.numeric(colnames(M))
    b  <- which.min(worst)
    cat(sprintf("  minimax-regret xi on %-4s = %.2f  (worst-case regret %.4g, %d cells%s)\n",
                lab, xi[b], worst[b], nrow(M),
                if (anyNA(worst)) sprintf(", %d xi dropped", sum(is.na(worst))) else ""))
    xi[b]
  }

  OPT <- list()      # per-cell RMSE optima, feeding the charge rule below
  MINIMAX <- list()  # per-arm minimax-regret xi, the fallback fixed value
  for (est in unique(D$estimator))
    for (ls in unique(D$lambda_src[D$estimator == est])) {
      sub0 <- D[D$estimator == est & D$lambda_src == ls & D$selector == "ebic", ]
      if (!nrow(sub0)) next
      cat(sprintf("=== %s  (lambda: %s) ===\n", est, ls))
      if (SHOW_CELLS)
        cat("    n  sigma  J |  RMSE-opt xi  [1 SE]  |  MCC-opt xi  [1 SE]  | at bound\n")
      RGr <- RGm <- list()
      for (cid in unique(cellid(sub0))) {
        s <- sub0[cellid(sub0) == cid, ]
        if (length(unique(s$rep)) < 10L) next
        a <- one_se(s, "rmse", FALSE); if (is.null(a)) next
        b <- if (s$J[1] > 0) one_se(s, "mcc", TRUE) else NULL
        RGr[[cid]] <- a$regret
        if (!is.null(b)) RGm[[cid]] <- b$regret        # J = 0 excluded
        OPT[[length(OPT) + 1L]] <- data.frame(
          estimator = est, lambda_src = ls, n = s$n[1], sigma = s$sigma[1],
          J = s$J[1], xi_rmse = a$best,
          xi_mcc = if (is.null(b)) NA_real_ else b$best,
          stringsAsFactors = FALSE)
        if (SHOW_CELLS)
          cat(sprintf("  %4d  %.2f  %d |  %.2f  [%.2f, %.2f]  |  %s |  %4.1f%%\n",
                      s$n[1], s$sigma[1], s$J[1], a$best, a$lo, a$hi,
                      if (is.null(b)) "   n/a (J=0)   " else
                        sprintf("%.2f  [%.2f, %.2f]", b$best, b$lo, b$hi),
                      100*mean(s$at_boundary)))
      }
      MINIMAX[[paste(est, ls)]] <- list(rmse = minimax(RGr, "RMSE"),
                                        mcc  = minimax(RGm, "MCC"))
      cat("\n")
    }
  OPT <- do.call(rbind, OPT)

  ## ---- false discoveries where there is nothing to find ------------------
  z0 <- D[D$J == 0 & D$selector == "ebic", ]
  if (nrow(z0)) {
    cat("=== false discoveries in the J = 0 cells ===\n")
    cat("  mean number of jumps declared when the truth has none\n")
    shown <- intersect(XI_GRID, c(0, 0.1, 0.2, 0.25, 0.3, 0.5, 1))
    ag <- aggregate(k_decl ~ estimator + lambda_src + xi + n,
                    z0[z0$xi %in% shown, ], function(v) mean(v, na.rm = TRUE))
    w <- reshape(ag, idvar = c("estimator", "lambda_src", "n"),
                 timevar = "xi", direction = "wide")
    print(w[order(w$estimator, w$lambda_src, w$n), ], row.names = FALSE, digits = 3)
    cat("\n")
  }

  ## ======================================================================
  ##  4. df per active warp, and the total charge per declared jump
  ## ======================================================================
  cat("=== effective degrees of freedom against active warps ===\n")
  slopes <- c()
  for (est in unique(FK$estimator)) {
    s <- FK[FK$estimator == est & is.finite(FK$df) & is.finite(FK$k_pkg), ]
    if (nrow(s) < 50L) next
    fit <- lm(df ~ k_pkg, data = s)
    slopes[est] <- unname(coef(fit)[2])
    cat(sprintf("  %-4s slope %6.3f (se %.3f): the plain BIC already levies\n",
                est, coef(fit)[2], summary(fit)$coefficients[2, 2]))
    cat(sprintf("       %.3f*log(n) per declared jump, before the EBIC term.\n",
                coef(fit)[2]))
  }
  cat("\n")

  ## The criterion charges  s*log(n)  through the plain BIC (df grows with the
  ## number of active warps) plus  2*xi*log(n-1)  through the EBIC term.  If
  ## the RMSE-optimal TOTAL charge is flat in n, then xi must fall with n and
  ## the rule  xi(n) = (C - s*log n) / (2 log(n-1))  follows.  Tested here
  ## rather than assumed: a significant trend in log(n) means the total charge
  ## is not constant and a fixed xi is the better description.
  cat("=== is the optimal total charge per jump constant in n? ===\n")
  RULE <- list()
  for (est in unique(OPT$estimator))
    for (ls in unique(OPT$lambda_src[OPT$estimator == est])) {
      o <- OPT[OPT$estimator == est & OPT$lambda_src == ls & OPT$J > 0 &
                 is.finite(OPT$xi_rmse), ]
      if (nrow(o) < 6L || is.na(slopes[est])) next
      s <- slopes[[est]]
      o$charge <- s*log(o$n) + 2*o$xi_rmse*log(o$n - 1)
      ft <- lm(charge ~ log(n), data = o)
      pv <- summary(ft)$coefficients[2, 4]
      ci <- tryCatch(confint(ft)[2, ], error = function(e) c(NA_real_, NA_real_))
      C  <- mean(o$charge)
      ## Declaring the charge constant on p > 0.05 alone would be reading
      ## absence of evidence as evidence of absence, and the test is weak when
      ## few cells feed it.  Require a reasonable number of cells and a CI that
      ## actually excludes a materially sloped charge; otherwise say so.
      const <- pv > 0.05 && nrow(o) >= 12L && is.finite(ci[1]) &&
               diff(ci) < 2 && abs(coef(ft)[2]) < 0.5
      RULE[[paste(est, ls)]] <- list(C = C, s = s, constant = const)
      cat(sprintf("  %-4s (lambda %s): mean charge %.2f, trend %+.3f per log(n) [%.2f, %.2f], p = %.3f, %d cells\n",
                  est, ls, C, coef(ft)[2], ci[1], ci[2], pv, nrow(o)))
      cat(if (const)
            "       flat enough -> prefer the xi(n) rule in the last column\n"
          else if (pv <= 0.05)
            "       significant trend -> the charge is NOT constant; use a fixed xi\n"
          else
            "       not resolved (wide interval or too few cells) -> use a fixed xi\n")
      byn <- aggregate(cbind(xi_rmse, charge) ~ n, o, mean)
      byn$xi_rule <- (C - s*log(byn$n)) / (2*log(byn$n - 1))
      print(byn, row.names = FALSE, digits = 3)
      cat("\n")
    }

  ## ======================================================================
  ##  5. Cross-validation against EBIC, at a fixed xi and at xi(n)
  ##
  ##  The fixed-xi arm is what the package does today.  Comparing CV only
  ##  against that conflates two questions, because the calibration above
  ##  shows the optimal xi moves with n: a constant xi is handicapped at large
  ##  n.  Both arms are therefore reported.
  ## ======================================================================
  A <- D[D$estimator == "atw", ]
  if (nrow(A) && length(RULE)) {
    r  <- RULE[["atw cv"]];    if (is.null(r))  r  <- RULE[["atw ml"]]
    mx <- MINIMAX[["atw cv"]]; if (is.null(mx)) mx <- MINIMAX[["atw ml"]]
    ## Use the n-dependent rule only where section 4 actually supports it;
    ## otherwise compare CV against the best fixed xi, which is the honest
    ## alternative rather than a rule the data did not license.
    use_rule <- !is.null(r) && isTRUE(r$constant)
    xi_fix   <- if (!is.null(mx) && is.finite(mx$rmse)) mx$rmse else XI_REF
    xi_of <- function(nn) if (!use_rule) xi_fix else
      XI_GRID[which.min(abs(XI_GRID - (r$C - r$s*log(nn))/(2*log(nn - 1))))]

    rows <- list()
    for (cid in unique(cellid(A))) {
      z <- A[cellid(A) == cid, ]
      xin <- xi_of(z$n[1])
      grab <- function(sel, ls, xi = NA) {
        q <- z[z$selector == sel & z$lambda_src == ls, ]
        if (!is.na(xi)) q <- q[abs(q$xi - xi) < 1e-9, ]
        if (!nrow(q)) return(NULL)
        data.frame(n = z$n[1], sigma = z$sigma[1], J = z$J[1],
                   lab = if (is.na(xi)) paste0(sel, " [", ls, "]")
                         else sprintf("ebic xi=%.2f [%s]", xi, ls),
                   rmse = mean(q$rmse, trim = ROBUST_TRIM, na.rm = TRUE),
                   rmse_med = median(q$rmse, na.rm = TRUE),
                   mcc = mean(q$mcc, na.rm = TRUE),
                   k = mean(q$k_decl, na.rm = TRUE), stringsAsFactors = FALSE)
      }
      for (ls in unique(z$lambda_src)) {
        rows <- c(rows, list(grab("ebic", ls, XI_REF)))
        if (abs(xin - XI_REF) > 1e-9) rows <- c(rows, list(grab("ebic", ls, xin)))
        for (sel in c("cv_M", "cv_lambda_M", "cv_lambda_ebic_M", "oracle"))
          rows <- c(rows, list(grab(sel, ls)))
      }
    }
    CV <- do.call(rbind, rows)
    fn <- file.path(OUT_DIR, "cv_comparison.csv")
    utils::write.csv(CV, fn, row.names = FALSE)

    cat("=== ssj_atw: CV against EBIC, at a fixed xi and at xi(n) ===\n")
    cat(sprintf("  second EBIC arm uses %s:  %s\n",
                if (use_rule) "the xi(n) rule from section 4"
                else "the minimax-regret fixed xi (section 4 did not license a rule)",
                paste(sprintf("n%d: %.2f", N_GRID, vapply(N_GRID, xi_of, 0)),
                      collapse = "   ")))
    cat(sprintf("  full 48-cell table written to %s\n", basename(fn)))
    cat("\n  averaged over sigma, cells with at least one jump (J >= 1):\n")
    s1 <- aggregate(cbind(rmse, mcc, k) ~ lab + n, CV[CV$J > 0, ], mean)
    print(s1[order(s1$n, s1$lab), ], row.names = FALSE, digits = 4)
    cat("\n  cells with no jumps (J = 0); k is the false-discovery count:\n")
    s0 <- aggregate(cbind(rmse, k) ~ lab + n, CV[CV$J == 0, ], mean)
    print(s0[order(s0$n, s0$lab), ], row.names = FALSE, digits = 4)
    cat("\n")
  }

  cat("Chen and Chen (2008) floor for p = n-1, xi >= 1 - log(n)/(2 log(n-1)):\n")
  for (nn in N_GRID) cat(sprintf("  n = %3d : %.3f\n", nn, 1 - log(nn)/(2*log(nn - 1))))
  cat("\nAdd the implicit charge from section 4 to the calibrated xi before\n")
  cat("comparing against this floor: the plain BIC already pays part of it.\n")

} else stop("argument must be \"probe\", \"run\" or \"analyse\"")
