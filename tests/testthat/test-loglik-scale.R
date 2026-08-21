## ssj_mle(), ssj_kfs() and ssj_atw() divide y by its standard deviation before
## fitting.  Until jumps 2.0.0 they reported the log-likelihood of that rescaled
## copy, which is invariant to the scale of y -- so multiplying y by 1000 left
## the reported log-likelihood, and every information criterion built on it,
## unchanged.  The likelihood of the data as supplied must instead fall by the
## Jacobian of the map, one factor of the scale per observation:
##
##      loglik(c * y)  =  loglik(y) - nobs * log(c).
##
## That identity is asserted below.  It needs no knowledge of the internals,
## and it fails under the old behaviour, where the difference was zero.

mk <- function(seed, n = 50) {
  set.seed(seed)
  x <- sort(runif(n))
  list(x = x, y = sin(2*pi*x) + 0.5*(x > 0.45) + rnorm(n, sd = 0.1))
}

test_that("the reported log-likelihood is on the scale of the data supplied", {
  d <- mk(9001); n <- length(d$y)
  for (cc in c(2, 10, 1000)) {
    ## lambda is a variance ratio, hence free of the scale of y: the same
    ## value is the right one for both data sets
    lam <- unname(ssj_mle(d$y, d$x, maxsum = 0)$pars["lambda"])

    k1 <- ssj_kfs(d$y,      d$x, lambda = lam, maxsum = 0)
    k2 <- ssj_kfs(cc*d$y,   d$x, lambda = lam, maxsum = 0)
    expect_equal(k2$loglik, k1$loglik - n*log(cc), tolerance = 1e-6)

    a1 <- ssj_atw(d$y,      d$x, lambda = lam, maxsum = 0)
    a2 <- ssj_atw(cc*d$y,   d$x, lambda = lam, maxsum = 0)
    expect_equal(a2$loglik, a1$loglik - n*log(cc), tolerance = 1e-6)
  }
})

test_that("ssj_mle's log-likelihood carries the same Jacobian", {
  d <- mk(9002); n <- length(d$y)
  m1 <- ssj_mle(d$y,      d$x, maxsum = 0)
  m2 <- ssj_mle(10*d$y,   d$x, maxsum = 0)
  expect_equal(m2$loglik, m1$loglik - n*log(10), tolerance = 1e-4)
  ## sigma is in units of y per x^{3/2}, so it scales with y; lambda is a
  ## variance ratio and does not
  expect_equal(unname(m2$pars["sigma"]), 10*unname(m1$pars["sigma"]),
               tolerance = 1e-3)
  expect_equal(unname(m2$pars["lambda"]), unname(m1$pars["lambda"]),
               tolerance = 1e-3)
})

test_that("the information criteria inherit the correction consistently", {
  d <- mk(9003); n <- length(d$y)
  lam <- unname(ssj_mle(d$y, d$x, maxsum = 0)$pars["lambda"])
  for (f in list(function(y) ssj_kfs(y, d$x, lambda = lam, maxsum = 0.2),
                 function(y) ssj_atw(y, d$x, lambda = lam, maxsum = 0.2))) {
    g1 <- f(d$y); g2 <- f(10*d$y)
    ## aic = 2*(df - loglik), so it must rise by 2*n*log(10)
    for (nm in c("aic", "bic", "hq")) {
      expect_equal(unname(g2$ic[nm]), unname(g1$ic[nm]) + 2*n*log(10),
                   tolerance = 1e-4)
    }
    expect_equal(g2$df, g1$df, tolerance = 1e-4)   # df is scale free
  }
})

test_that("rescaling y does not change which model is selected", {
  ## the Jacobian is the same at every maxsum, so it cancels in the comparison
  ## the auto_* functions make: this is why the correction is safe
  d <- mk(9004)
  lam <- unname(ssj_mle(d$y, d$x, maxsum = 0)$pars["lambda"])
  g <- seq(0, 1, length.out = 6)
  a1 <- auto_ssj_atw(d$y,    d$x, lambda = lam, grid = g, ic = "ebic")
  a2 <- auto_ssj_atw(50*d$y, d$x, lambda = lam, grid = g, ic = "ebic")
  expect_equal(a1$maxsum, a2$maxsum)
  k1 <- auto_ssj_kfs(d$y,    d$x, lambda = lam, grid = g, ic = "ebic")
  k2 <- auto_ssj_kfs(50*d$y, d$x, lambda = lam, grid = g, ic = "ebic")
  expect_equal(k1$maxsum, k2$maxsum)
})
