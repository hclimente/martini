X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))

test_that("output is as expected", {
  
  cones <- sigmod_(X, Y, test_gwas$map$snp.name, test_gi, .47, .0183)
  cones.cv <- sigmod.cv_(X, Y, test_gwas$map$snp.name, test_gi)
  
  expect_true(all(V(cones) == V(cones.cv)))
  expect_true(all(E(cones) == E(cones.cv)))
  
})
