test_that("output is as expected", {
  
  cones <- scones(test_gwas, test_gi, 10, 1)
  cones.cv <- scones.cv(test_gwas, test_gi)
  
  expect_true(all(V(cones) == V(cones.cv)))
  expect_true(all(E(cones) == E(cones.cv)))
  
})
