test_that("output is as expected", {
  
  cones <- sigmod(test_gwas, test_gi, 7.37, 0.202)
  cones.cv <- sigmod.cv(test_gwas, test_gi)
  
  expect_true(all(V(cones) == V(cones.cv)))
  expect_true(all(E(cones) == E(cones.cv)))
  
})
