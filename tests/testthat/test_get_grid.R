library(martini)

set.seed(0)

test_that("default args are being set", {
  expect_equal(get_grid(c = 1)$etas, rep(1, 10))
  expect_equal(get_grid(c = 1)$lambdas, 
              10^seq(-1, 1, length.out = 10), tolerance=0.01)
})

test_that("default we can change values", {
  expect_equal(get_grid(c = 1, etas = c(3,4,5))$etas, c(3,4,5))
  expect_error(get_grid(etas = c("a","b"))$etas)
  expect_equal(get_grid(c = 1, lambdas = c(3,4,5))$lambdas, c(3,4,5))
  expect_error(get_grid(etas = c(1,2,3), lambdas = c("a","b")),
               'specify a valid lambdas or an association vector.')
  expect_error(get_grid(etas = c("a","b"), lambdas = c(1,2,3)),
               'specify a valid etas or an association vector.')
  expect_error(get_grid(c = 1, debug = 3))
})
