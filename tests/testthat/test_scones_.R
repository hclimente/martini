X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))

test_that("output is as expected", {
  
  cones <- scones_(X, Y, test_gwas$map$snp.name, test_gi, 10, 1)
  cones.cv <- scones.cv_(X, Y, test_gwas$map$snp.name, test_gi)
  
  expect_true(all(V(cones) == V(cones.cv)))
  expect_true(all(E(cones) == E(cones.cv)))
  
})
