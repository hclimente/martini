library(martini)
data(examplegwas)

set.seed(0)

test_that("the ordering of the map is not relevant", {
  expect_equal(test_snp_clusters(map, net, N = 100), test_snp_clusters(map[sample(nrow(map)),], net, N = 100))
})