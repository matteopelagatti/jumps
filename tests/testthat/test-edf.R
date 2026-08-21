## The effective degrees of freedom reported by the Kalman-based estimators is
## sum(w), where w is the diagonal of the hat matrix computed inside llt() and
## llt_delta() with the Koopman and Harvey (2003) weight recursion.  With the
## variance parameters held fixed the map y -> smoothed_level is exactly linear,
## so that diagonal can be obtained by central differences and compared term by
## term.  These tests pin the recursion down: before jumps 2.0.0 the index on N
## was off by one and the weights could fall outside (0, 1] and sum to a
## negative number.

## Run the smoother at fixed parameters and return both the fitted signal and
## the weights the C++ routine reports.
smooth_at <- function(y, tau, sig2, lambda, wt = rep(1, length(y))) {
  n <- length(y); n1 <- n + 1L
  var_eps <- lambda * sig2 / wt
  a1 <- rep(y[!is.na(y)][1L], n1); a2 <- numeric(n1)
  p11 <- rep(1e5, n1); p12 <- numeric(n1); p22 <- rep(1e5, n1)
  k1 <- numeric(n); k2 <- numeric(n); inn <- numeric(n); fv <- numeric(n)
  r1 <- numeric(n1); r2 <- numeric(n1)
  n11 <- numeric(n1); n12 <- numeric(n1); n22 <- numeric(n1)
  e <- numeric(n1); d <- numeric(n1); w <- numeric(n)
  llt_delta(y, tau, var_eps, sig2 * tau^3 / 3, sig2 * tau, sig2 * tau^2 / 2,
            a1, a2, p11, p12, p22, k1, k2, inn, fv,
            r1, r2, n11, n12, n22, e, d, w)
  list(mu = (a1 + p11 * r1 + p12 * r2)[1:n], w = w)
}

## Exact hat-matrix diagonal by central differences.
hat_diag <- function(y, tau, sig2, lambda, wt = rep(1, length(y)), h = 1e-5) {
  vapply(seq_along(y), function(j) {
    if (is.na(y[j])) return(0)
    yp <- y; yp[j] <- yp[j] + h
    ym <- y; ym[j] <- ym[j] - h
    ((smooth_at(yp, tau, sig2, lambda, wt)$mu -
      smooth_at(ym, tau, sig2, lambda, wt)$mu) / (2 * h))[j]
  }, numeric(1))
}

set.seed(11)
n  <- 25
x  <- sort(runif(n))
y  <- sin(2 * pi * x) + 0.6 * (x > 0.5) + rnorm(n, sd = 0.15)
tau <- c(diff(x), 1)

test_that("the reported weights are the hat-matrix diagonal", {
  for (lambda in c(0.01, 0.4, 5)) {
    got <- smooth_at(y, tau, 0.09, lambda)$w
    expect_equal(got, hat_diag(y, tau, 0.09, lambda), tolerance = 1e-5)
  }
})

test_that("the weight recursion is correct with observation weights", {
  wt  <- c(3, 1, 5, 2, 1, 1, 4, 2, 6, 1, 3, 2, 1, 2, 1, 4, 2, 1, 3, 1, 2, 5, 1, 2, 3)
  got <- smooth_at(y, tau, 0.09, 0.4, wt)$w
  expect_equal(got, hat_diag(y, tau, 0.09, 0.4, wt), tolerance = 1e-5)
})

test_that("a missing observation gets zero weight", {
  ym <- y; ym[c(1, 7, 8, n)] <- NA          # including the first
  got <- smooth_at(ym, tau, 0.09, 0.4)$w
  expect_equal(got[c(1, 7, 8, n)], rep(0, 4))
  expect_equal(got, hat_diag(ym, tau, 0.09, 0.4), tolerance = 1e-5)
})

test_that("weights lie in (0, 1] and df is between 2 and n", {
  fits <- list(
    ssj_mle(faithful$waiting, faithful$eruptions, maxsum = 1),
    ssj_kfs(faithful$waiting, faithful$eruptions, lambda = 1, maxsum = 1),
    ssj_fix(faithful$waiting, faithful$eruptions, lambda = 1, maxsum = 2),
    hpfj(as.numeric(Nile), maxsum = 50),
    hpfj_fix(as.numeric(Nile), lambda = 1600, maxsum = 50)
  )
  for (f in fits) {
    expect_gte(f$df, 2 - 1e-6)
    expect_lte(f$df, f$nobs + 1e-6)
  }
  hf <- hpfj(as.numeric(Nile), maxsum = 50)
  expect_true(all(hf$weights > 0))
  expect_true(all(hf$weights <= 1 + 1e-10))
  expect_equal(sum(hf$weights), hf$df)
})

test_that("llt() and llt_delta() agree on an equally spaced design", {
  n <- 20
  yy <- as.numeric(Nile)[1:n]
  n1 <- n + 1L
  arg <- function() list(a1 = rep(yy[1], n1), a2 = numeric(n1),
                         p11 = rep(1e5, n1), p12 = numeric(n1), p22 = rep(1e5, n1),
                         k1 = numeric(n), k2 = numeric(n), inn = numeric(n),
                         fv = numeric(n), r1 = numeric(n1), r2 = numeric(n1),
                         n11 = numeric(n1), n12 = numeric(n1), n22 = numeric(n1),
                         e = numeric(n1), d = numeric(n1), w = numeric(n))
  ve <- rep(1600, n); vn <- rep(1/3, n); vz <- rep(1, n); vc <- rep(1/2, n)
  a <- arg()
  llt(yy, ve, vn, vz, vc, a$a1, a$a2, a$p11, a$p12, a$p22, a$k1, a$k2,
      a$inn, a$fv, a$r1, a$r2, a$n11, a$n12, a$n22, a$e, a$d, a$w)
  b <- arg()
  llt_delta(yy, rep(1, n), ve, vn, vz, vc, b$a1, b$a2, b$p11, b$p12, b$p22,
            b$k1, b$k2, b$inn, b$fv, b$r1, b$r2, b$n11, b$n12, b$n22,
            b$e, b$d, b$w)
  expect_equal(a$w, b$w)
})
