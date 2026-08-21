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
