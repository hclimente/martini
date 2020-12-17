library(igraph)
library(martini)

test_that("output is as expected", {
  
  cones <- sigmod.cv(test_gwas, test_gi)
  
  expect_equal(dim(cones), dim(test_gwas$map) + c(0,3))
  expect_equal(class(cones), "data.frame")
  expect_equal(class(cones$selected), "logical")
  expect_equal(class(cones$c), "numeric")
  expect_equal(class(cones$module), "numeric")
  
})

test_that("we recover causal SNPs", {
  
  cones <- sigmod.cv(test_gwas, test_gi, etas = 0, lambdas = 0)
  
  # wrong eta and lambda return the trivial solution
  expect_equal(sum(cones$selected), sum(cones$c > 0))
  
  set.seed(42)
  cones <- sigmod.cv(test_gwas, test_gi, 
                     etas = seq(2, 0, length=10), 
                     lambdas = seq(2, 0, length=10))
  
  skip_on_os("windows")
  expect_equal(sum(cones$selected), sum(grepl("[AC]", cones$snp)))
  
  scores <- martini:::snp_test(test_gwas, data.frame(), 'chi2')
  c <- cones$c
  names(c) <- cones$snp
  expect_equal(c, scores)
  
  set.seed(42)
  cones <- sigmod.cv(test_gwas, test_gi, 
                     etas = seq(2, 0, length=10), 
                     lambdas = seq(2, 0, length=10),
                     criterion = 'bic')
  
  skip_on_os("windows")
  expect_equal(sum(cones$selected), sum(grepl("[AC]", cones$snp)))
  
  c <- cones$c
  names(c) <- cones$snp
  expect_equal(c[cones$selected], scores[grepl("[AC]", names(scores))])
  
})
