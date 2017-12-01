library(igraph)
library(martini)

source("big_network.R")

test_that("output is as expected", {
  
  cones <- search_cones(gwas, gi)
  
  expect_equal(dim(cones), dim(gwas$map) + c(0,3))
  expect_equal(class(cones), "data.frame")
  expect_equal(class(cones$selected), "logical")
  expect_equal(class(cones$c), "numeric")
  expect_equal(class(cones$module), "numeric")
  
})

test_that("we recover causal SNPs", {
  
  cones <- search_cones(gwas, gi, etas = 1, lambdas = 2)
  
  expect_equal(sum(cones$selected), nrow(cones))
  
  cones <- search_cones(gwas, gi)
  
  expect_equal(sum(cones$selected[sol]), pCausal)
  expect_equal(cones$c[sol], rep(96.15385, pCausal), tolerance = 1e-5)
  
})