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
  expect_equal(fit$nobs, length(faithful$waiting))
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

  fc <- ssj_cssd_potts(y, x, p = 0.999, gamma = 4)
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

  out_cssd <- capture.output(print(ssj_cssd_potts(faithful$waiting,
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
  f <- fit_atw_()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  before <- par("mfrow")
  expect_null(plot(f))
  expect_identical(par("mfrow"), before)
  expect_null(plot(f, show_warp = FALSE))
  expect_null(plot(f, thresh = 0.1))
  expect_null(plot(f, main = "a title of my own"))
  expect_null(plot(ssj_mle(faithful$waiting[1:60], faithful$eruptions[1:60],
                           maxsum = 1), prob = 0.95))
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
  expect_null(plot(f))                       # untitled
  expect_null(plot(f, main = "custom"))      # titled
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
