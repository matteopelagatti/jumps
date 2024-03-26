test_that("trigseas works", {
  expect_equal(
    trigseas(12, 4),
    cbind(cos1 = cos((0:11)*2*pi/4),
          cos2 = cos((0:11)*4*pi/4),
          sin1 = sin((0:11)*2*pi/4))
  )
})

test_that("dummyseas works", {
          expect_equal(
            dummyseas(12, 4),
            cbind(D1 = rep(c(1, -1/3, -1/3, -1/3), 3),
                  D2 = rep(c(-1/3, 1, -1/3, -1/3), 3),
                  D3 = rep(c(-1/3, -1/3, 1, -1/3), 3)
            )
          )
})

test_that("mse works", {
  expect_equal(mse(1:3, 2:4), 1)
})