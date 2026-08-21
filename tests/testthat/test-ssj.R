test_that("ssj_atw returns fitted values aligned with its own x slot", {
  fit <- ssj_atw(y = faithful$waiting, x = faithful$eruptions,
                 lambda = 1, maxsum = 150)

  expect_length(fit$smoothed_level, length(fit$x))
  expect_equal(fit$nobs, length(unique(faithful$eruptions)))
  expect_true(all(is.finite(fit$smoothed_level)))
  expect_length(fit$gamma, fit$nobs - 1L)
})

test_that("ssj_atw refuses missing values in y", {
  x <- sort(runif(30))
  y <- sin(4*pi*x)
  y[5] <- NA
  expect_error(ssj_atw(y, x, lambda = 1, maxsum = 0.5), "missing")
})

test_that("ssj_kfs and ssj_mle tolerate missing values in y", {
  set.seed(3)
  x <- sort(runif(40))
  y <- sin(2*pi*x) + rnorm(40, sd = 0.1)
  y[7] <- NA

  fit_kfs <- ssj_kfs(y, x, lambda = 1, maxsum = 0.3)
  expect_true(all(is.finite(fit_kfs$smoothed_level)))
  expect_length(fit_kfs$smoothed_level, length(x))

  fit_mle <- ssj_mle(y, x, maxsum = 0.3)
  expect_true(all(is.finite(fit_mle$smoothed_level)))
  expect_true(fit_mle$pars[["lambda"]] > 0)
})

test_that("ssj_fix returns the documented slots", {
  fit <- ssj_fix(y = faithful$waiting, x = faithful$eruptions,
                 lambda = 1, maxsum = 150)

  expect_true(all(c("opt", "nobs", "df", "maxsum", "loglik", "pars", "sigmas",
                    "ic", "smoothed_level", "var_smoothed_level", "x") %in%
                    names(fit)))
  # faithful$eruptions has ties, which are now collapsed before fitting
  expect_equal(fit$nobs, length(unique(faithful$eruptions)))
  expect_length(fit$x, fit$nobs)
  expect_true(all(fit$var_smoothed_level >= 0))
})

# ---- the ssj class and its methods ---------------------------------------

fit_atw_ <- function()
  ssj_atw(faithful$waiting, faithful$eruptions, lambda = 1, maxsum = 1)

test_that("the ssj_* fitters return classed objects without losing slots", {
  fa <- faithful[seq(1, nrow(faithful), by = 4), ]
  fits <- list(
    atw = fit_atw_(),
    mle = ssj_mle(fa$waiting, fa$eruptions, maxsum = 1),
    kfs = ssj_kfs(fa$waiting, fa$eruptions, lambda = 1, maxsum = 1),
    fix = ssj_fix(fa$waiting, fa$eruptions, lambda = 1, maxsum = 150)
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_s3_class(f, "ssj")
    expect_true(is.list(f))
    expect_identical(f$method, nm)
    expect_true(all(c("nobs", "df", "loglik", "ic", "smoothed_level", "x") %in%
                      names(f)))
    expect_length(f$y, length(f$x))
  }
})

test_that("auto_* wrappers keep the class and record their own call", {
  f <- auto_ssj_atw(faithful$waiting, faithful$eruptions, lambda = 1,
                    grid = c(0, 0.5, 1))
  expect_s3_class(f, "ssj")
  expect_match(paste(deparse(f$call), collapse = ""), "auto_ssj_atw")
})

test_that("ssj_breaks dates a break at the first observation of the new regime", {
  f <- fit_atw_()
  b <- ssj_breaks(f)
  expect_s3_class(b, "data.frame")
  expect_named(b, c("index", "x", "warp", "jump"))
  skip_if(nrow(b) == 0L)

  # the largest warp sits on (x[t], x[t+1]) and must be dated at t+1
  t_max <- which.max(f$gamma)
  expect_identical(as.integer(b$index[which.max(b$warp)]), as.integer(t_max + 1L))
  expect_equal(b$x, f$x[b$index])
  expect_equal(b$jump, f$smoothed_level[b$index] - f$smoothed_level[b$index - 1L])
})

test_that("the detection rule requires both of the paper's conditions", {
  f <- fit_atw_()
  g <- f$gamma
  b <- ssj_breaks(f, frac_total = 0.05, frac_delta = 0.30)
  keep <- g > 0.05*sum(g) & g > 0.30*mean(diff(f$x))
  expect_equal(nrow(b), sum(keep))

  # an absolute threshold overrides both conditions
  expect_equal(nrow(ssj_breaks(f, thresh = 0)), sum(g > 0))
  expect_equal(nrow(ssj_breaks(f, thresh = max(g))), 0L)
})

test_that("ssj_breaks handles the CSSD competitor and TV denoising", {
  set.seed(1)
  x <- sort(runif(80)); y <- sin(4*pi*x) + (x > 0.5) + rnorm(80, sd = 0.1)

  fc <- ssj_cssd_dp(y, x, p = 0.999, gamma = 4)
  expect_s3_class(fc, "ssj")
  expect_identical(fc$method, "cssd")
  bc <- ssj_breaks(fc)
  expect_equal(nrow(bc), fc$n_disc)
  expect_true(all(is.na(bc$warp)))

  ft <- ssj_cssd(y, x, lambda = 2)
  expect_identical(ft$method, "cssd_tv")
  # the ADMM iterate is only approximately sparse: the share-of-total rule must
  # keep the break set small rather than flagging nearly every interval
  expect_lt(nrow(ssj_breaks(ft)), 15L)
})

test_that("print.ssj does not leak partial matches and labels the budget", {
  fa <- faithful[seq(1, nrow(faithful), by = 4), ]
  out_mle <- capture.output(print(ssj_mle(fa$waiting, fa$eruptions, maxsum = 1)))
  # `x$p` used to partial-match `pars` and print a bogus "p" parameter
  expect_false(any(grepl("^  p +=", out_mle)))
  expect_true(any(grepl("sigma", out_mle)))
  expect_true(any(grepl("domain warps", out_mle)))

  out_fix <- capture.output(
    print(ssj_fix(fa$waiting, fa$eruptions, lambda = 1, maxsum = 150)))
  expect_true(any(grepl("additional standard deviations", out_fix)))

  out_cssd <- capture.output(print(ssj_cssd_dp(faithful$waiting,
                                                 faithful$eruptions,
                                                 p = 0.99, gamma = 5)))
  expect_true(any(grepl("Storath and Weinmann", out_cssd)))
})

test_that("print.ssj returns its argument invisibly", {
  f <- fit_atw_()
  expect_output(res <- withVisible(print(f)))
  expect_false(res$visible)
  expect_identical(res$value, f)
})

test_that("plot.ssj draws and restores the graphical parameters", {
  # this is about the base-graphics path: it splits the device into two panels
  # and must put par() back as it found it
  f <- fit_atw_()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  before <- par("mfrow")
  expect_null(plot(f, use_ggplot = FALSE))
  expect_identical(par("mfrow"), before)
  expect_null(plot(f, use_ggplot = FALSE, show_warp = FALSE))
  expect_null(plot(f, use_ggplot = FALSE, thresh = 0.1))
  expect_null(plot(f, use_ggplot = FALSE, main = "a title of my own"))
  expect_null(plot(ssj_mle(faithful$waiting[1:60], faithful$eruptions[1:60],
                           maxsum = 1), prob = 0.95, use_ggplot = FALSE))
  expect_identical(par("mfrow"), before)
})

test_that("ssj_breaks rejects objects of the wrong class", {
  expect_error(ssj_breaks(list(a = 1)), "ssj")
})

test_that("the default maxsum is the range of x, not the rescaled one", {
  # The default refers to x, and R forces it lazily.  ssj_kfs()/ssj_mle()
  # rescale x to [0, 1] internally, so without an explicit force() the promise
  # would evaluate to 1 rather than to the range of x.
  fa  <- faithful[seq(1, nrow(faithful), by = 4), ]
  rng <- max(fa$eruptions) - min(fa$eruptions)
  expect_gt(rng, 1)                       # otherwise the test proves nothing

  expect_equal(ssj_mle(fa$waiting, fa$eruptions)$maxsum, rng)
  expect_equal(ssj_kfs(fa$waiting, fa$eruptions, lambda = 1)$maxsum, rng)
  expect_equal(ssj_atw(fa$waiting, fa$eruptions, lambda = 1)$maxsum, rng)

  # the three warping estimators share one default, equal to the upper end of
  # the grid searched by the corresponding auto_* wrapper
  expect_identical(formals(ssj_mle)$maxsum, formals(ssj_kfs)$maxsum)
  expect_identical(formals(ssj_mle)$maxsum, formals(ssj_atw)$maxsum)

  # and the constraint is actually enforced at that value
  f <- ssj_mle(fa$waiting, fa$eruptions)
  expect_lte(sum(f$gamma), f$maxsum * (1 + 1e-8))

  # an explicit budget equal to the default reproduces the default fit
  expect_equal(ssj_mle(fa$waiting, fa$eruptions, maxsum = rng)$smoothed_level,
               f$smoothed_level)

  # the default must not depend on the order in which x is supplied
  set.seed(11); i <- sample(nrow(fa))
  expect_equal(ssj_mle(fa$waiting[i], fa$eruptions[i])$maxsum, rng)
})

test_that("plot.ssj draws no title unless one is asked for", {
  # the method labels are long, so the default must be no title at all
  expect_null(eval(formals(jumps:::plot.ssj)$main))
  f <- fit_atw_()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  # base graphics: drawn, so the value is NULL either way
  expect_null(plot(f, use_ggplot = FALSE))                  # untitled
  expect_null(plot(f, use_ggplot = FALSE, main = "custom"))  # titled
  # ggplot2 is the default and returns the object, so the title can be read off
  skip_if_not_installed("ggplot2")
  expect_null(ggplot2::ggplot_build(plot(f))$plot$labels$title)
  expect_identical(ggplot2::ggplot_build(plot(f, main = "custom"))$plot$labels$title,
                   "custom")
})

test_that("fits do not depend on the order in which observations are supplied", {
  # var(diff(y)) is not invariant to row order and scales the diffuse prior on
  # the slope, so it must be computed after sorting by x.  With distinct x the
  # sorted (x, y) are identical under any input permutation, hence so must the
  # fit be.
  set.seed(4)
  x <- sort(runif(50)); y <- sin(3*pi*x) + rnorm(50, sd = 0.1)
  i <- sample(length(x))
  expect_false(isTRUE(all.equal(var(diff(y)), var(diff(y[i])))))  # premise

  expect_equal(ssj_mle(y[i], x[i], maxsum = 0.3)$smoothed_level,
               ssj_mle(y,    x,    maxsum = 0.3)$smoothed_level)
  expect_equal(ssj_kfs(y[i], x[i], lambda = 1, maxsum = 0.3)$smoothed_level,
               ssj_kfs(y,    x,    lambda = 1, maxsum = 0.3)$smoothed_level)
  expect_equal(ssj_fix(y[i], x[i], lambda = 1, maxsum = 0.3)$smoothed_level,
               ssj_fix(y,    x,    lambda = 1, maxsum = 0.3)$smoothed_level)
  expect_equal(ssj_atw(y[i], x[i], lambda = 1, maxsum = 0.3)$smoothed_level,
               ssj_atw(y,    x,    lambda = 1, maxsum = 0.3)$smoothed_level)
})

test_that("predict.ssj returns the cubic the model implies", {
  set.seed(2026)
  x <- sort(runif(50)); y <- sin(4*pi*x) + (x > 0.5) + rnorm(50, sd = 0.1)
  fit <- ssj_kfs(y, x, lambda = 0.01, maxsum = 0.3)

  # at the sampling sites it must reproduce the fit exactly
  expect_equal(predict(fit, fit$x), fit$smoothed_level)
  expect_equal(predict(fit), fit$smoothed_level)

  # order of the query points is preserved
  xq <- c(0.8, 0.2, 0.5, 0.35)
  expect_equal(predict(fit, xq)[order(xq)], predict(fit, sort(xq)))

  # the cubic is NOT the linear interpolation of the fitted values
  xg  <- seq(min(x), max(x), length.out = 300)
  cub <- predict(fit, xg)
  lin <- stats::approx(fit$x, fit$smoothed_level, xout = xg, rule = 2)$y
  expect_gt(max(abs(cub - lin)), 1e-6)
  expect_true(all(is.finite(cub)))

  # beyond the data the extension is linear: second differences vanish
  xo <- seq(max(x) + 0.01, max(x) + 0.20, length.out = 12)
  expect_lt(max(abs(diff(diff(predict(fit, xo))))), 1e-8)
})

test_that("predict.ssj works for every estimator", {
  set.seed(5)
  x <- sort(runif(40)); y <- sin(3*pi*x) + rnorm(40, sd = 0.1)
  fits <- list(
    mle  = ssj_mle(y, x, maxsum = 0.2),
    kfs  = ssj_kfs(y, x, lambda = 1, maxsum = 0.2),
    atw  = ssj_atw(y, x, lambda = 1, maxsum = 0.2),
    fix  = ssj_fix(y, x, lambda = 1, maxsum = 5),
    cssd = ssj_cssd_dp(y, x, p = 0.99, gamma = 2),
    tv   = ssj_cssd(y, x, lambda = 0.5)
  )
  for (nm in names(fits)) {
    f <- fits[[nm]]
    expect_equal(predict(f, f$x), f$smoothed_level, tolerance = 1e-6,
                 info = nm)
    expect_true(all(is.finite(predict(f, seq(0, 1, length.out = 50)))), info = nm)
  }
})

test_that("CSSD extrapolates linearly beyond a segment, not with the end cubic", {
  set.seed(1)
  x <- sort(runif(60)); y <- sin(4*pi*x) + (x > 0.5) + rnorm(60, sd = 0.1)
  f <- ssj_cssd_dp(y, x, p = 0.999, gamma = 4)

  # outside the data range the second difference of the prediction must vanish
  for (xo in list(seq(min(x) - 0.30, min(x) - 0.02, length.out = 10),
                  seq(max(x) + 0.02, max(x) + 0.30, length.out = 10))) {
    p <- predict(f, xo)
    expect_lt(max(abs(diff(diff(p)))), 1e-8)
  }
})

test_that("extrapolation is not amplified by a warp on a boundary interval", {
  # The warps live strictly BETWEEN observations, so outside the sample the
  # dilation is 1.  Extrapolating with the physical slope just inside the end
  # interval would multiply the slope by tau/delta, which on a boundary
  # interval carrying a warp can be a factor of tens.
  set.seed(2026)
  x <- sort(runif(60)); y <- sin(4*pi*x) + (x > 0.5) + rnorm(60, sd = 0.1)
  f <- ssj_kfs(y, x, lambda = 0.01, maxsum = 0.3)
  n <- length(f$x)

  h <- 1e-5
  left  <- (predict(f, f$x[1] - h) - predict(f, f$x[1] - 2*h)) / h
  right <- (predict(f, f$x[n] + 2*h) - predict(f, f$x[n] + h)) / h
  expect_equal(left,  f$smoothed_slope[1], tolerance = 1e-5)
  expect_equal(right, f$smoothed_slope[n], tolerance = 1e-5)

  # and it stays bounded: extrapolating a tenth of the range must not move the
  # prediction by more than a few times the range of the fit
  rng <- diff(range(f$smoothed_level))
  far <- predict(f, c(f$x[1] - 0.1*diff(range(f$x)),
                      f$x[n] + 0.1*diff(range(f$x))))
  expect_true(all(abs(far - mean(f$smoothed_level)) < 20 * rng))
})

test_that("weighted collapsing reproduces the fit on the raw replicates", {
  # Collapsing ties to their mean and passing the replicate counts as weights
  # is sufficient at a FIXED smoothing parameter: the discarded within-group
  # sum of squares does not involve f.
  f <- function(x) sin(2*pi*x) + 0.6*(x > 0.5)
  set.seed(11)
  xu <- sort(runif(30)); m <- sample(1:5, 30, replace = TRUE)
  x  <- rep(xu, times = m); y <- f(x) + rnorm(length(x), sd = 0.15)
  ybar <- as.numeric(tapply(y, x, mean)); w <- as.numeric(table(x))
  expect_gt(max(m), min(m))            # weights must actually be unequal

  # the fitters collapse ties themselves, so both fits live on the 40 sites
  raw <- ssj_kfs(y, x, lambda = 0.05, maxsum = 0)
  wgt <- ssj_kfs(ybar, xu, lambda = 0.05, maxsum = 0, weights = w)
  expect_length(raw$smoothed_level, length(xu))
  expect_equal(raw$smoothed_level, wgt$smoothed_level, tolerance = 1e-3)

  raw <- ssj_atw(y, x, lambda = 0.05, maxsum = 0)
  wgt <- ssj_atw(ybar, xu, lambda = 0.05, maxsum = 0, weights = w)
  expect_equal(raw$smoothed_level, wgt$smoothed_level)
  expect_equal(raw$df, wgt$df)
})

test_that("unit weights reproduce the unweighted fit exactly", {
  set.seed(5)
  x <- sort(runif(40)); y <- sin(3*pi*x) + rnorm(40, sd = 0.1)
  o <- rep(1, 40)
  for (nm in c("kfs", "mle", "fix", "atw")) {
    a <- switch(nm, kfs = ssj_kfs(y, x, lambda = 1, maxsum = 0.2),
                    mle = ssj_mle(y, x, maxsum = 0.2),
                    fix = ssj_fix(y, x, lambda = 1, maxsum = 2),
                    atw = ssj_atw(y, x, lambda = 1, maxsum = 0.2))
    b <- switch(nm, kfs = ssj_kfs(y, x, lambda = 1, maxsum = 0.2, weights = o),
                    mle = ssj_mle(y, x, maxsum = 0.2, weights = o),
                    fix = ssj_fix(y, x, lambda = 1, maxsum = 2, weights = o),
                    atw = ssj_atw(y, x, lambda = 1, maxsum = 0.2, weights = o))
    expect_equal(a$smoothed_level, b$smoothed_level, info = nm)
    expect_equal(a$ic, b$ic, info = nm)
  }
})

test_that("weights are validated and forwarded by the auto_* wrappers", {
  set.seed(6)
  x <- sort(runif(30)); y <- sin(3*pi*x) + rnorm(30, sd = 0.1)
  expect_error(ssj_kfs(y, x, lambda = 1, weights = c(-1, rep(1, 29))), "positive")
  expect_error(ssj_kfs(y, x, lambda = 1, weights = rep(1, 5)), "same length")
  expect_silent(f <- ssj_kfs(y, x, lambda = 1, maxsum = 0.2, weights = 3))
  expect_s3_class(f, "ssj")

  w <- rep(c(1, 4), 15)
  for (nm in c("atw", "kfs")) {
    a <- switch(nm, atw = auto_ssj_atw(y, x, lambda = 1, grid = c(0, 0.2)),
                    kfs = auto_ssj_kfs(y, x, lambda = 1, grid = c(0, 0.2)))
    b <- switch(nm, atw = auto_ssj_atw(y, x, lambda = 1, grid = c(0, 0.2), weights = w),
                    kfs = auto_ssj_kfs(y, x, lambda = 1, grid = c(0, 0.2), weights = w))
    expect_false(isTRUE(all.equal(a$smoothed_level, b$smoothed_level)), info = nm)
  }
})

test_that("ss_within restores the noise information lost by collapsing", {
  f <- function(x) sin(2*pi*x) + 0.6*(x > 0.5)
  set.seed(11)
  xu <- sort(runif(40)); m <- sample(1:6, 40, replace = TRUE)
  x  <- rep(xu, times = m); y <- f(x) + rnorm(length(x), sd = 0.15)
  yb <- as.numeric(tapply(y, x, mean)); w <- as.numeric(table(x))
  ssw <- sum(tapply(y, x, function(v) sum((v - mean(v))^2)))
  ve <- function(f) f$pars[["lambda"]] * f$pars[["sigma"]]^2

  raw  <- ssj_mle(y,  x,  maxsum = 0)
  bare <- ssj_mle(yb, xu, maxsum = 0, weights = w)
  full <- ssj_mle(yb, xu, maxsum = 0, weights = w, ss_within = ssw)

  # collapsing alone overstates the noise; ss_within recovers the raw value
  expect_gt(ve(bare) / ve(raw), 2)
  expect_equal(ve(full), ve(raw), tolerance = 1e-3)
  expect_equal(full$pars[["lambda"]], raw$pars[["lambda"]], tolerance = 1e-2)

  # and the pure-smooth fit then matches the raw one
  expect_equal(full$smoothed_level, raw$smoothed_level, tolerance = 1e-3)
})

test_that("ss_within is validated and inert when omitted", {
  set.seed(7)
  x <- sort(runif(30)); y <- sin(3*pi*x) + rnorm(30, sd = 0.1)
  a <- ssj_mle(y, x, maxsum = 0.2)
  b <- ssj_mle(y, x, maxsum = 0.2, ss_within = NULL)
  expect_equal(a$smoothed_level, b$smoothed_level)

  # needs weights summing to more than length(y), i.e. genuine replicate counts
  expect_error(ssj_mle(y, x, maxsum = 0.2, ss_within = 1), "weights")
  expect_error(ssj_mle(y, x, maxsum = 0.2, weights = rep(2, 30),
                       ss_within = -1), "non-negative")
})

test_that("tied x are collapsed automatically by every estimator", {
  f <- function(x) sin(2*pi*x) + 0.6*(x > 0.5)
  set.seed(11)
  xu <- sort(runif(40)); m <- sample(1:6, 40, replace = TRUE)
  x  <- rep(xu, times = m); y <- f(x) + rnorm(length(x), sd = 0.15)
  yb <- as.numeric(tapply(y, x, mean)); w <- as.numeric(table(x))
  ssw <- sum(tapply(y, x, function(v) sum((v - mean(v))^2)))

  for (nm in c("mle", "kfs", "fix", "atw")) {
    auto <- switch(nm, mle = ssj_mle(y, x, maxsum = 0),
                       kfs = ssj_kfs(y, x, lambda = 0.05, maxsum = 0),
                       fix = ssj_fix(y, x, lambda = 0.05, maxsum = 0),
                       atw = ssj_atw(y, x, lambda = 0.05, maxsum = 0))
    expect_equal(length(auto$x), length(xu), info = nm)
    expect_equal(auto$x, xu, tolerance = 1e-8, info = nm)

    man <- switch(nm,
      mle = ssj_mle(yb, xu, maxsum = 0, weights = w, ss_within = ssw),
      kfs = ssj_kfs(yb, xu, lambda = 0.05, maxsum = 0, weights = w),
      fix = ssj_fix(yb, xu, lambda = 0.05, maxsum = 0, weights = w),
      atw = ssj_atw(yb, xu, lambda = 0.05, maxsum = 0, weights = w))
    # collapsing internally must equal collapsing by hand and passing weights
    expect_equal(auto$smoothed_level, man$smoothed_level,
                 tolerance = 1e-3, info = nm)
  }
})

test_that("a missing first observation no longer poisons the filter", {
  # llt()/llt_delta() initialised t = 0 without the missing-value guard used in
  # the loop, so a leading NA made the whole log-likelihood NA.
  set.seed(3)
  x <- sort(runif(50)); y <- sin(3*pi*x) + rnorm(50, sd = 0.1)
  y1 <- y; y1[1] <- NA
  for (fit in list(ssj_mle(y1, x, maxsum = 0),
                   ssj_kfs(y1, x, lambda = 1, maxsum = 0),
                   ssj_fix(y1, x, lambda = 1, maxsum = 2))) {
    expect_true(all(is.finite(fit$smoothed_level)))
    expect_true(is.finite(fit$loglik))
  }
  y2 <- y; y2[c(1, 2, 25, 50)] <- NA
  expect_true(all(is.finite(ssj_mle(y2, x, maxsum = 0.2)$smoothed_level)))
})

test_that("a wrong-length parinit is rejected with a message naming the length", {
  # A wrong-length parinit used to reach nloptr as an NA-padded x0 and fail with
  # the opaque "x0 contains NA".
  set.seed(5)
  x <- sort(runif(30)); y <- sin(2*pi*x) + rnorm(30, sd = 0.1)

  # ssj_mle expects n+1 = 31: (sigma, lambda, gamma_1, ..., gamma_{n-1})
  expect_error(ssj_mle(y, x, maxsum = 0.2, parinit = c(0.1, 1, rep(0, 10))),
               "length 31, not 12")
  expect_error(ssj_mle(y, x, maxsum = 0.2, parinit = c(0.1, 1, rep(0, 10))),
               "gamma_1")

  # ssj_kfs expects n = 30: (sigma, gamma_1, ..., gamma_{n-1})
  expect_error(ssj_kfs(y, x, lambda = 1, maxsum = 0.2, parinit = rep(0.1, 31)),
               "length 30, not 31")

  # ssj_fix expects n+1 = 31: (sigma, sigma_1, ..., sigma_n)
  expect_error(ssj_fix(y, x, lambda = 1, maxsum = 2, parinit = rep(0.1, 29)),
               "length 31, not 29")

  # a correct-length parinit still goes through
  expect_silent(ssj_kfs(y, x, lambda = 1, maxsum = 0.2,
                        parinit = c(sd(y)/10, rep(0, 29))))
})

test_that("the parinit length error counts distinct x, and says so", {
  # The expected length follows the number of DISTINCT x after ties have been
  # collapsed, so length(y) is the wrong thing for a caller to reason from.
  x <- rep(sort(runif(20)), each = 3)          # 60 observations, 20 distinct x
  set.seed(7)
  y <- sin(2*pi*x) + rnorm(60, sd = 0.1)

  # length(y) + 1 = 61 looks right but is not: the fitters want 21 / 20
  expect_error(ssj_mle(y, x, maxsum = 0.2, parinit = rep(0.1, 61)),
               "length 21, not 61")
  expect_error(ssj_mle(y, x, maxsum = 0.2, parinit = rep(0.1, 61)),
               "distinct values of 'x'")
  expect_error(ssj_kfs(y, x, lambda = 1, maxsum = 0.2, parinit = rep(0.1, 60)),
               "length 20, not 60")
  expect_error(ssj_fix(y, x, lambda = 1, maxsum = 2, parinit = rep(0.1, 61)),
               "length 21, not 61")
})
