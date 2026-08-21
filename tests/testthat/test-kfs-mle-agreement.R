## ssj_kfs() fits exactly the model ssj_mle() fits, with lambda supplied rather
## than estimated.  Handed ssj_mle()'s own lambda, it must therefore reach the
## same optimum.  Before jumps 2.0.0 it did not: its single fixed starting
## value for sigma (sdy/10) could sit orders of magnitude from the optimum,
## and the optimizer diverged in the great majority of cases -- returning
## sigma ~ 1e6 and a log-likelihood hundreds of units below the attainable one.
## ssj_fix() shared the defect. These tests pin the agreement down.

mk <- function(seed, n = 60) {
  set.seed(seed)
  x <- sort(runif(n))
  list(x = x,
       y = sin(2*pi*x) + 0.5*sin(5*pi*x) + 0.5*(x > 0.4) - 0.4*(x > 0.7) +
           rnorm(n, sd = 0.08))
}

test_that("ssj_kfs at the MLE lambda is never worse than ssj_mle", {
  ## This is the regression guard.  The failure it protects against was
  ## ssj_kfs landing hundreds of log-likelihood units BELOW ssj_mle on the
  ## same data at the same lambda.  Equality is not asserted, because the
  ## likelihood in sigma is multimodal for some data sets and either fitter
  ## can settle on the better mode -- at seed 4002 it is ssj_mle that reaches
  ## the worse one.  What must never happen again is ssj_kfs falling short.
  for (s in 1:6) {
    d <- mk(4000 + s)
    m <- ssj_mle(d$y, d$x, maxsum = 0)
    k <- ssj_kfs(d$y, d$x, lambda = unname(m$pars["lambda"]), maxsum = 0)
    expect_gte(k$loglik, m$loglik - 1e-4)
  }
})

test_that("where both reach the same optimum, the fits coincide", {
  for (s in 1:6) {
    d <- mk(4000 + s)
    m <- ssj_mle(d$y, d$x, maxsum = 0)
    k <- ssj_kfs(d$y, d$x, lambda = unname(m$pars["lambda"]), maxsum = 0)
    if (abs(k$loglik - m$loglik) > 1e-5) next     # different modes: skip
    expect_equal(unname(k$pars["sigma"]), unname(m$pars["sigma"]), tolerance = 1e-3)
    expect_equal(k$df, m$df, tolerance = 1e-4)
    expect_equal(k$smoothed_level, m$smoothed_level, tolerance = 1e-4)
  }
})

test_that("the agreement survives a positive warping budget", {
  d <- mk(4010)
  m <- ssj_mle(d$y, d$x, maxsum = 0.2)
  k <- ssj_kfs(d$y, d$x, lambda = unname(m$pars["lambda"]), maxsum = 0.2)
  expect_equal(k$loglik, m$loglik, tolerance = 1e-3)
  expect_equal(sum(k$gamma), sum(m$gamma), tolerance = 1e-3)
  expect_equal(k$smoothed_level, m$smoothed_level, tolerance = 1e-3)
})

test_that("the agreement is invariant to the scale of x and y", {
  d <- mk(4020)
  for (fx in c(1, 10)) for (fy in c(1, 5)) {
    m <- ssj_mle(d$y*fy, d$x*fx, maxsum = 0)
    k <- ssj_kfs(d$y*fy, d$x*fx, lambda = unname(m$pars["lambda"]), maxsum = 0)
    expect_gte(k$loglik, m$loglik - 1e-4)
    expect_equal(k$smoothed_level, m$smoothed_level, tolerance = 1e-4)
  }
})

test_that("ssj_fix at maxsum = 0 is the same model as ssj_kfs", {
  ## With every additional standard deviation forced to zero, ssj_fix has
  ## var_eta = sigma^2 delta^3/3, var_zeta = sigma^2 delta,
  ## cov = sigma^2 delta^2/2 and var_eps = lambda sigma^2 / w: exactly what
  ## ssj_kfs fits at gamma = 0.  The fitted curve, the scale parameter and the
  ## degrees of freedom must therefore agree.
  ##
  ## The reported log-likelihoods do NOT agree, and deliberately are not
  ## compared here: ssj_mle and ssj_kfs rescale y internally and report the
  ## likelihood of the rescaled copy, so they sit below the likelihood of the
  ## data as supplied by the Jacobian n*log(sd(y)) (about 12.5 units on this
  ## example), whereas ssj_fix does not rescale.  The offset is constant in
  ## maxsum, so it cancels inside auto_ssj_*(); see the note in NEWS.md.
  for (s in 1:3) {
    d <- mk(4030 + s)
    m <- ssj_mle(d$y, d$x, maxsum = 0)
    lam <- unname(m$pars["lambda"])
    k <- ssj_kfs(d$y, d$x, lambda = lam, maxsum = 0)
    f <- ssj_fix(d$y, d$x, lambda = lam, maxsum = 0)
    expect_equal(f$smoothed_level, k$smoothed_level, tolerance = 1e-4)
    expect_equal(f$df, k$df, tolerance = 1e-4)
    expect_equal(unname(f$pars["sigma"]), unname(k$pars["sigma"]), tolerance = 0.05)
  }
})

test_that("neither fitter returns an absurd sigma", {
  ## the failure mode was a runaway scale parameter, so bound it against the
  ## data: sigma is a slope-disturbance standard deviation per unit of x
  for (s in 1:3) {
    d <- mk(4040 + s)
    m <- ssj_mle(d$y, d$x, maxsum = 0)
    lam <- unname(m$pars["lambda"])
    for (f in list(ssj_kfs(d$y, d$x, lambda = lam, maxsum = 0),
                   ssj_fix(d$y, d$x, lambda = lam, maxsum = 0))) {
      expect_true(is.finite(f$pars["sigma"]))
      expect_lt(unname(f$pars["sigma"]), 1e4)
      expect_gt(unname(f$pars["sigma"]), 1e-8)
    }
  }
})
