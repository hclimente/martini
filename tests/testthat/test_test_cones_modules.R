load("examplegwas.rda")
set.seed(0)

cones <- examplegwas$cones
net <- examplegwas$net

test_that("the ordering of the map is not relevant", {
  
  shuffled <- cones[sample(nrow(cones)),]
  
  expect_equal(test_cones_modules(cones, net, nperm = 100), 
               test_cones_modules(shuffled, net, nperm = 100))

})
