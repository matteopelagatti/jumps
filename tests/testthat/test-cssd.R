make_step <- function(n = 60, seed = 1) {
  set.seed(seed)
  x <- sort(runif(n))
  y <- sin(4 * pi * x) + (x > 0.5) + rnorm(n, sd = 0.05)
  list(x = x, y = y)
}

test_that("ssj_cssd_dp recovers a single discontinuity", {
  d   <- make_step()
  fit <- ssj_cssd_dp(d$y, d$x, p = 0.999, gamma = 4)

  expect_length(fit$smoothed_level, length(d$x))
  expect_true(all(is.finite(fit$smoothed_level)))
  expect_equal(fit$n_disc, length(fit$disc_locs))
  expect_gte(fit$n_disc, 1L)
  expect_true(any(abs(fit$disc_locs - 0.5) < 0.05))
})

test_that("ssj_cssd_dp validates its arguments", {
  d <- make_step(20)
  expect_error(ssj_cssd_dp(d$y, d$x[-1]), "same length")
  expect_error(ssj_cssd_dp(d$y, d$x, p = 1), "\\(0, 1\\)")
  expect_error(ssj_cssd_dp(d$y, d$x, gamma = 0), "positive")
  expect_error(ssj_cssd_dp(replace(d$y, 1, NA), d$x), "NAs")
})

test_that("predict_cssd_dp preserves the order of xnew", {
  d    <- make_step()
  fit  <- ssj_cssd_dp(d$y, d$x, p = 0.999, gamma = 4)
  xq   <- c(0.9, 0.1, 0.5, 0.3)

  p_shuffled <- predict_cssd_dp(fit, xq)
  p_sorted   <- predict_cssd_dp(fit, sort(xq))

  expect_equal(p_shuffled[order(xq)], p_sorted)
  # at the data sites the prediction must reproduce the fit
  expect_equal(predict_cssd_dp(fit, fit$x), fit$smoothed_level,
               tolerance = 1e-6)
})

test_that("predict_cssd_dp covers the gap around a discontinuity", {
  d   <- make_step()
  fit <- ssj_cssd_dp(d$y, d$x, p = 0.999, gamma = 4)
  skip_if(fit$n_disc == 0L)

  dl   <- fit$disc_locs[1]
  # a dense grid straddling the discontinuity: no query point may be dropped
  xq   <- seq(dl - 0.02, dl + 0.02, length.out = 25)
  pred <- predict_cssd_dp(fit, xq)

  expect_true(all(is.finite(pred)))
  expect_true(all(abs(pred) > .Machine$double.eps))
  expect_true(max(abs(diff(pred))) > 0.5)   # the jump is visible
})

test_that("auto_ssj_cssd works with the default (regular) x", {
  d   <- make_step(40)
  set.seed(2)
  fit <- auto_ssj_cssd(d$y)

  expect_equal(fit$x, seq_along(d$y))
  expect_named(fit$cv_par, c("p", "gamma"))
  expect_true(fit$cv_par[["p"]] > 0 && fit$cv_par[["p"]] < 1)
  expect_true(fit$cv_par[["gamma"]] > 0)
})

test_that("ssj_cssd validates lambda and detects the jump", {
  d <- make_step()
  expect_error(ssj_cssd(d$y, d$x, lambda = -1), "positive")
  expect_error(ssj_cssd(d$y, d$x, lambda = c(1, 2)), "positive")

  fit <- ssj_cssd(d$y, d$x, lambda = 2)
  expect_length(fit$gamma, fit$nobs - 1L)
  expect_true(all(is.finite(fit$smoothed_level)))
})

test_that("tied x values are collapsed consistently", {
  d <- make_step(40)
  x <- c(d$x, d$x[10])
  y <- c(d$y, d$y[10])
  fit <- ssj_cssd_dp(y, x, p = 0.999, gamma = 4)
  expect_equal(fit$nobs, 40L)
  expect_length(fit$smoothed_level, 40L)
})

## --- the tuning search must not depend on the units of y --------------------
## gamma is a charge per discontinuity set against the data-fit and smoothness
## terms, both of which scale with the square of the units of y.  auto_ssj_cssd()
## therefore searches on a standardised response and maps gamma back; before it
## did so, the same data in MWh and in GWh gave 146 discontinuities and none.

test_that("auto_ssj_cssd is invariant to the units of y", {
  set.seed(4242)
  x <- sort(runif(60))
  y <- sin(2*pi*x) + 0.8*(x > 0.55) + rnorm(60, sd = 0.15)

  ref <- NULL
  for (cc in c(1, 1e3, 1e6)) {
    set.seed(1)
    f <- suppressWarnings(auto_ssj_cssd(y = y/cc, x = x))
    if (is.null(ref)) {
      ref <- list(lev = f$smoothed_level, nd = f$n_disc,
                  p = f$cv_par[["p"]], gz = f$cv_par[["gamma"]] / stats::var(y))
      next
    }
    expect_equal(f$smoothed_level * cc, ref$lev, tolerance = 1e-6, info = cc)
    expect_equal(f$n_disc, ref$nd, info = cc)
    expect_equal(f$cv_par[["p"]], ref$p, tolerance = 1e-6, info = cc)
    # gamma must move with var(y): the ratio is the scale-free quantity
    expect_equal(f$cv_par[["gamma"]] / stats::var(y/cc), ref$gz,
                 tolerance = 1e-6, info = cc)
  }
})

test_that("auto_ssj_cssd reports gamma in the units of y", {
  set.seed(4243)
  x <- sort(runif(50)); y <- sin(2*pi*x) + rnorm(50, sd = 0.2)
  set.seed(1)
  f <- suppressWarnings(auto_ssj_cssd(y = 1000*y, x = x))
  set.seed(1)
  g <- suppressWarnings(auto_ssj_cssd(y = y, x = x))
  # a 1000-fold change in y is a 1e6-fold change in gamma
  expect_equal(f$cv_par[["gamma"]] / g$cv_par[["gamma"]], 1e6,
               tolerance = 1e-4)
  # and cv_score, being a squared error, carries the same factor
  expect_equal(f$cv_score / g$cv_score, 1e6, tolerance = 1e-3)
})

test_that("a jump penalty at the edge of the search range warns", {
  # A grid deliberately confined far below the optimum: the refinement escapes
  # it, and the selected penalty then sits outside the range that was searched,
  # which is exactly the situation the user has to be told about -- it is how
  # the electricity-demand fit came to report 146 discontinuities in silence.
  set.seed(4244)
  x <- sort(runif(50)); y <- sin(2*pi*x) + rnorm(50, sd = 0.2)
  expect_warning(auto_ssj_cssd(y = y, x = x, gamma_grid = c(1e-6, 2e-6)),
                 "edge of its search range")
  # the default grid, on the same data, has no such complaint
  expect_no_warning(auto_ssj_cssd(y = y, x = x))
})

test_that("gamma_grid is validated", {
  set.seed(4245)
  x <- sort(runif(30)); y <- sin(2*pi*x) + rnorm(30, sd = 0.2)
  expect_error(auto_ssj_cssd(y = y, x = x, gamma_grid = c(-1, 1)), "positive")
  expect_error(auto_ssj_cssd(y = y, x = x, gamma_grid = numeric(0)), "positive")
})

test_that("an extreme gamma_grid does not overflow the q parameterisation", {
  # gamma = p*q/(1-q): reaching gamma ~ 1e13 needs logit(q) ~ 35, past which
  # 1 - q underflows and gamma becomes Inf.  A grid in that region must still
  # return a usable fit (the classical spline) rather than fail.
  set.seed(4246)
  x <- sort(runif(40)); y <- sin(2*pi*x) + rnorm(40, sd = 0.2)
  f <- suppressWarnings(auto_ssj_cssd(y = y, x = x, gamma_grid = c(1e12, 1e15)))
  expect_s3_class(f, "ssj")
  expect_true(all(is.finite(f$smoothed_level)))
  expect_true(is.finite(f$cv_par[["gamma"]]))
  expect_equal(f$n_disc, 0L)          # such a penalty forbids discontinuities
})
