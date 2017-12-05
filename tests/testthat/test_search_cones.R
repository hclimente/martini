library(igraph)
library(martini)

test_that("output is as expected", {
  
  cones <- search_cones(gwas, gi)
  
  expect_equal(dim(cones), dim(minigwas$map) + c(0,3))
  expect_equal(class(cones), "data.frame")
  expect_equal(class(cones$selected), "logical")
  expect_equal(class(cones$c), "numeric")
  expect_equal(class(cones$module), "numeric")
  
})

test_that("we recover causal SNPs", {
  
  gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
  cones <- search_cones(minigwas, gi, etas = 1, lambdas = 2)
  
  # wrong eta and lambda return the trivial solution
  expect_equal(sum(cones$selected), nrow(cones))
  
  cones <- search_cones(minigwas, gi)
  
  expect_equal(sum(cones$selected), sum(grepl("[AC]", cones$snp)))
  expect_equal(cones$c[cones$selected], rep(96.15385, sum(grepl("[AC]", cones$snp))), tolerance = 1e-5)
  
})