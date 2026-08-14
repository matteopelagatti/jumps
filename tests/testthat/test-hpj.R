sim_hp <- function(n = 60, seed = 202311) {
  set.seed(seed)
  mu <- 100*cos(3*pi/n*(1:n)) - ((1:n) > n/2)*n
  mu + rnorm(n, sd = 20)
}

test_that("hpfj_fix reports gamma from the right parameter", {
  y  <- sim_hp()
  hf <- hpfj_fix(y, lambda = 60, maxsum = 60)

  # parameter vector is (sd_zeta, sqrt_gamma, sigma_eta_1, ..., sigma_eta_n)
  expect_length(hf$opt$solution, length(y) + 2L)
  expect_equal(unname(hf$pars["gamma"]), hf$opt$solution[2]^2)
  expect_equal(unname(hf$pars["sigma_noise"]),
               hf$opt$solution[1]*sqrt(60))
})

test_that("BIC.hpj accepts objects with extra classes", {
  y  <- sim_hp()
  h1 <- hpj(y, 50)
  h2 <- hpj(y, 60)
  h3 <- h2
  class(h3) <- c("hpj", "extra")

  expect_type(BIC(h1), "double")
  tab <- BIC(h1, h2, h3)
  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 3L)
  expect_named(tab, c("df", "BIC"))
})

test_that("hpj rejects xreg and invalid scl", {
  y <- sim_hp(30)
  expect_error(hpj(y, 50, xreg = matrix(1, 30, 1)), "xreg")
  expect_error(hpj(y, 50, scl = -1), "scl")
})

test_that("logLik and nobs methods are consistent", {
  y  <- sim_hp()
  hp <- hpj(y, 50)
  ll <- logLik(hp)

  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "nobs"), nobs(hp))
  expect_equal(nobs(hp), length(y))
})

test_that("the plot method returns a base plot invisibly", {
  y  <- sim_hp(30)
  hp <- hpj(y, 50)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_null(plot(hp, use_ggplot = FALSE))
  expect_null(plot(hp, prob = 0.95, use_ggplot = FALSE))
})
