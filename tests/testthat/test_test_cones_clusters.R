library(martini)
data(examplegwas)

set.seed(0)

test_that("the ordering of the map is not relevant", {
  expect_equal(test_cones_clusters(examplegwas$map, examplegwas$net, N = 100), 
               test_cones_clusters(examplegwas$map[sample(nrow(examplegwas$map)),], examplegwas$net, N = 100))
})