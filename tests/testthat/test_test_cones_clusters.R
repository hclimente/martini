library(martini)
data(examplegwas)

set.seed(0)

test_that("the ordering of the map is not relevant", {
  expect_equal(test_cones_clusters(examplegwas$cones, examplegwas$net, N = 100), 
               test_cones_clusters(examplegwas$cones[sample(nrow(examplegwas$cones)),], examplegwas$net, N = 100))
})