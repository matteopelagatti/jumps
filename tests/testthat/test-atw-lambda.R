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

## --- tied covariate values -------------------------------------------------
## The fold partition must group replicates: ssj_atw() collapses tied x inside
## whatever sample it is handed, so an observation-level split leaves siblings
## of a held-out point in the training set and the criterion goes optimistic.

mk_tied <- function(seed, nsite = 24, reps = 2:6) {
  set.seed(seed)
  xu <- round(seq(0, 1, length.out = nsite), 4)
  # NB: not sample(reps, ...) -- a length-1 'reps' would make sample() draw
  # from seq_len(reps) instead of from reps itself
  m  <- if (length(reps) == 1L) rep(as.integer(reps), nsite)
        else sample(reps, nsite, replace = TRUE)
  x  <- rep(xu, times = m)
  list(x = x, y = sin(2*pi*x) + 0.6*(x > 0.45) + rnorm(length(x), sd = 0.1),
       site = rep(seq_along(xu), times = m), nsite = nsite)
}

test_that("default folds never split the replicates of a tied x", {
  d <- mk_tied(8101)
  expect_gt(length(d$x), d$nsite)             # ties must actually be present

  # the partition is internal, so probe it through the seed: with grouping, a
  # fit is reproducible from a seed and the recorded cv surface is finite
  for (s in 1:5) {
    set.seed(s)
    f <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT,
                                          grid = G6, nfolds = 3))
    expect_true(any(is.finite(f$cv_surface)))
  }

  # and the grouped partition this now builds is exactly what passing
  # site-level folds by hand produces: same seed, same answer
  set.seed(42)
  auto <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT,
                                           grid = G6, nfolds = 3))
  set.seed(42)
  sf <- sample(rep_len(seq_len(3L), d$nsite))
  byhand <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT,
                                             grid = G6, folds = sf[d$site]))
  expect_equal(auto$cv_surface, byhand$cv_surface)
  expect_equal(auto$pars[["lambda"]], byhand$pars[["lambda"]])
})

test_that("nfolds is bounded by the number of distinct x, not of observations", {
  d <- mk_tied(8102, nsite = 6, reps = 4L)   # 24 observations, 6 sites
  expect_equal(length(d$x), 24L)
  expect_silent(f <- suppressWarnings(
    auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6, nfolds = 6)))
  expect_error(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6,
                               nfolds = 7),
               "number of distinct values")
})

test_that("a supplied partition that splits a site warns", {
  d <- mk_tied(8103)
  bad <- sample(rep_len(1:3, length(d$x)))     # observation-level: splits sites
  expect_warning(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6,
                                 folds = bad),
                 "splits the replicates")
  good <- sample(rep_len(1:3, d$nsite))[d$site]
  expect_no_warning(suppressWarnings({   # boundary warning may still fire
    withCallingHandlers(
      auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT, grid = G6, folds = good),
      warning = function(w) if (grepl("splits the replicates", conditionMessage(w)))
        stop("site-splitting warning raised for a grouped partition"))
  }))
})

test_that("an untied design is unaffected by the grouping", {
  d <- mk(8104, n = 40)                        # continuous x: no ties
  expect_equal(length(unique(d$x)), length(d$x))
  set.seed(9)
  a <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT,
                                        grid = G6, nfolds = 5))
  set.seed(9)
  b <- suppressWarnings(auto_ssj_atw_cv(d$y, d$x, lambda_grid = SHORT,
                                        grid = G6, nfolds = 5))
  expect_equal(a$cv_surface, b$cv_surface)
  expect_false(any(is.na(a$cv_surface)))
})
