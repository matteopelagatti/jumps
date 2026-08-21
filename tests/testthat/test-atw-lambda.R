## auto_ssj_atw_gcv() and auto_ssj_atw_cv() choose lambda without recourse to
## ssj_mle(), then hand it to auto_ssj_atw() for the budget.

mk <- function(seed, n = 60) {
  set.seed(seed)
  x <- sort(runif(n))
  list(x = x, y = sin(2*pi*x) + 0.6*(x > 0.45) + rnorm(n, sd = 0.1))
}
SHORT <- 10^seq(-6, -3, length.out = 4)
G6    <- seq(0, 1, length.out = 6)

test_that("both routines return a usable ssj fit", {
  d <- mk(8001)
  ## SHORT is deliberately narrow, so a boundary warning is expected here and
  ## is not what these assertions are about
  g <- suppressWarnings(auto_ssj_atw_gcv(d$y, d$x, lambda_grid = SHORT, grid = G6))
  c1 <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6, nfolds = 3))
  for (f in list(g, c1)) {
    expect_s3_class(f, "ssj")
    expect_equal(length(f$smoothed_level), length(f$x))
    expect_true(all(is.finite(f$smoothed_level)))
    expect_true(is.finite(f$pars["lambda"]))
    expect_true(f$maxsum %in% G6)
  }
  expect_identical(g$lambda_source, "gcv")
  expect_identical(c1$lambda_source, "cv")
})

test_that("the selected lambda comes from the supplied grid", {
  d <- mk(8002)
  g <- suppressWarnings(auto_ssj_atw_gcv(d$y, d$x, lambda_grid = SHORT, grid = G6))
  expect_true(unname(g$pars["lambda"]) %in% SHORT)
  expect_equal(nrow(g$lambda_path), length(SHORT))
  expect_equal(unname(g$pars["lambda"]),
               g$lambda_path$lambda[which.min(g$lambda_path$gcv)])

  cc <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6, nfolds = 3))
  expect_true(unname(cc$pars["lambda"]) %in% SHORT)
  expect_equal(dim(cc$cv_surface), c(length(SHORT), length(G6)))
  expect_equal(unname(cc$pars["lambda"]),
               cc$lambda_path$lambda[which.min(cc$lambda_path$cv)])
})

test_that("the default lambda grid follows the scale of x", {
  ## lambda has the units of x^3, so the same data on a stretched x must give a
  ## grid, and a selected lambda, scaled by that factor cubed
  d <- mk(8003)
  a <- suppressWarnings(auto_ssj_atw_gcv(d$y,   d$x,    grid = G6))
  b <- suppressWarnings(auto_ssj_atw_gcv(d$y, 10*d$x, grid = 10*G6))
  expect_equal(unname(b$pars["lambda"]), unname(a$pars["lambda"]) * 1000,
               tolerance = 1e-6)
  expect_equal(b$smoothed_level, a$smoothed_level, tolerance = 1e-4)
})

test_that("supplied folds make the cross-validation reproducible", {
  d <- mk(8004)
  fd <- rep_len(1:3, length(d$y))
  a <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6, folds = fd))
  b <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6, folds = fd))
  expect_equal(a$pars["lambda"], b$pars["lambda"])
  expect_equal(a$cv_surface, b$cv_surface)
})

test_that("a boundary solution is reported rather than passed off silently", {
  d <- mk(8005)
  ## a grid whose upper end is far too smooth forces the choice onto it
  expect_warning(auto_ssj_atw_gcv(d$y, d$x, grid = G6,
                                  lambda_grid = 10^seq(-12, -10, length.out = 3)),
                 "largest lambda|smallest lambda")
})

test_that("arguments are validated", {
  d <- mk(8006)
  expect_error(auto_ssj_atw_gcv(d$y, d$x[-1]), "same length")
  expect_error(auto_ssj_atw_gcv(c(NA, d$y[-1]), d$x), "missing values")
  expect_error(auto_ssj_atw_gcv(d$y, d$x, lambda_grid = c(-1, 1e-4)), "positive")
  expect_error(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, nfolds = 1),
               "between 2 and")
  expect_error(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, folds = 1:3),
               "same length")
})

test_that("weights are accepted by both routines", {
  d <- mk(8007, n = 50)
  w <- rep(c(1, 2), length.out = length(d$y))
  g <- suppressWarnings(auto_ssj_atw_gcv(d$y, d$x, weights = w, lambda_grid = SHORT, grid = G6))
  c1 <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, weights = w, lambda_grid = SHORT,
                                         grid = G6, nfolds = 3))
  expect_true(all(is.finite(g$smoothed_level)))
  expect_true(all(is.finite(c1$smoothed_level)))
})
