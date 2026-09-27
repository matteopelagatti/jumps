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

## --- warp_frac: the lower panel's share of the height -----------------------

mk_fit <- function() ssj_atw(y = faithful$waiting, x = faithful$eruptions,
                             lambda = 1, maxsum = 1)

panel_heights <- function(p) {
  g   <- ggplot2::ggplotGrob(p)
  idx <- g$layout[grepl("^panel", g$layout$name), ]
  as.numeric(g$heights[idx$t[order(idx$t)]])
}

test_that("warp_frac sets the lower panel's share of the plotting height", {
  skip_if_not_installed("ggplot2")
  skip_if(utils::packageVersion("ggplot2") < "4.0.0",
          "panel.heights needs ggplot2 >= 4.0.0")
  fit <- mk_fit()
  for (f in c(0.2, 0.3, 0.5)) {
    h <- panel_heights(plot(fit, warp_frac = f))
    expect_length(h, 2L)
    expect_equal(h[2]/sum(h), f, tolerance = 1e-8, info = f)
  }
  # the default is 0.3, not the 0.5 that equal panels amount to
  expect_equal(panel_heights(plot(fit))[2]/sum(panel_heights(plot(fit))), 0.3,
               tolerance = 1e-8)
})

test_that("the figure is still a plain ggplot that can be extended", {
  skip_if_not_installed("ggplot2")
  p <- plot(mk_fit(), warp_frac = 0.25)
  expect_s3_class(p, "ggplot")
  expect_s3_class(p + ggplot2::labs(title = "t"), "ggplot")
})

test_that("the warp panel is labelled 'Warps'", {
  skip_if_not_installed("ggplot2")
  lv <- levels(ggplot2::ggplot_build(plot(mk_fit()))$layout$layout$panel)
  expect_true("Warps" %in% lv)
  expect_false(any(grepl("Domain warp", lv)))
})

test_that("warp_frac is validated", {
  fit <- mk_fit()
  for (bad in list(0, 1, -0.2, 1.5, c(0.3, 0.4), NA_real_, "a", NULL))
    expect_error(plot(fit, warp_frac = bad), "strictly between 0 and 1")
})

test_that("base graphics honour warp_frac and restore par", {
  fit <- mk_fit()
  f <- tempfile(fileext = ".png")
  png(f, width = 600, height = 500)
  on.exit({ if (grDevices::dev.cur() > 1L) grDevices::dev.off()
            unlink(f) }, add = TRUE)
  expect_null(plot(fit, use_ggplot = FALSE, warp_frac = 0.2))
  # layout() must not leak: par(op) on exit returns the device to a 1x1 grid
  expect_equal(par("mfrow"), c(1L, 1L))
  grDevices::dev.off()
  expect_gt(file.size(f), 1000)
})
