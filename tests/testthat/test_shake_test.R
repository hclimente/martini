library(martini)
data(examplegwas)

set.seed(0)

test_that("the ordering of the map is not relevant", {
  expect_equal(shake_test(map, net, N = 100), shake_test(map[sample(nrow(map)),], net, N = 100))
})