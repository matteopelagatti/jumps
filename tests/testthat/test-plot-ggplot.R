## plot.ssj() gains the use_ggplot switch that plot.hpj() already had, and
## defaults to it, so the method returns a ggplot object unless asked otherwise.

fit_atw <- ssj_atw(y = faithful$waiting, x = faithful$eruptions,
                   lambda = 1, maxsum = 1)

test_that("ggplot2 is the default and base graphics returns NULL invisibly", {
  expect_identical(eval(formals(jumps:::plot.ssj)$use_ggplot), TRUE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_null(plot(fit_atw, use_ggplot = FALSE))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(fit_atw), "ggplot")
})

test_that("use_ggplot returns a ggplot object that prints", {
  skip_if_not_installed("ggplot2")
  p <- plot(fit_atw, use_ggplot = TRUE)
  expect_s3_class(p, "ggplot")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(print(p))
  expect_s3_class(p + ggplot2::labs(title = "x"), "ggplot")
})

test_that("the ggplot path honours show_warp, show_breaks, prob and main", {
  skip_if_not_installed("ggplot2")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  expect_silent(print(plot(fit_atw, use_ggplot = TRUE, show_warp = FALSE)))
  expect_silent(print(plot(fit_atw, use_ggplot = TRUE, show_breaks = FALSE)))
  expect_silent(print(plot(fit_atw, use_ggplot = TRUE, main = "title")))
  ## a band needs pointwise variances, which ssj_atw does not return
  m <- ssj_mle(faithful$waiting[1:60], faithful$eruptions[1:60], maxsum = 1)
  expect_silent(print(plot(m, use_ggplot = TRUE, prob = 0.95)))
})

test_that("the two panels are facets of one object, not two plots", {
  skip_if_not_installed("ggplot2")
  p <- plot(fit_atw, use_ggplot = TRUE, show_warp = TRUE)
  expect_s3_class(p, "ggplot")
  b <- ggplot2::ggplot_build(p)
  expect_gt(length(unique(b$data[[1]]$PANEL)), 0L)
})
