library(martini)
load("examplegwas.rda")
set.seed(0)

test_that("the ordering of the map is not relevant", {
  expect_equal(test_cones_modules(examplegwas$cones, examplegwas$net, nperm = 100), 
               test_cones_modules(examplegwas$cones[sample(nrow(examplegwas$cones)),], examplegwas$net, nperm = 100))
})