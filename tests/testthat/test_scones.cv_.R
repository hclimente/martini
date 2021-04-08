X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))

test_that("output is as expected", {
  
  cones <- scones.cv_(X, Y, test_gwas$map$snp.name, test_gi)
  
  expect_equal(class(cones), "igraph")
  expect_true(all(names(V(cones)) %in% test_map$snp))
  
})

test_that("we recover causal SNPs", {
  
  cones <- scones.cv_(X, Y, test_gwas$map$snp.name, test_gi)
  
  skip_on_os("windows")
  expect_equal(gorder(cones), sum(grepl("[AC]", test_map$snp)))
  expect_equal(sort(names(V(cones))), 
               sort(grep("[AC]", test_map$snp, value = TRUE)))
  
})