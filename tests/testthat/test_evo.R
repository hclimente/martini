library(martini)

source("big_network.R")

test_that("output is as expected", {
  
  settings <- get_evo_settings(etas = 1, lambdas = 2)
  test <- evo(X, Y, W, settings)
  
  expect_equal(length(test), 4)
  expect_equal(class(test), "list")
  expect_equal(test$eta, 1)
  expect_equal(test$lambda, 2)
  expect_equal(class(test$selected), "numeric")
  expect_equal(class(test$c), "numeric")
  
})

test_that("we recover causal SNPs", {
  
  settings <- get_evo_settings(etas = 1, lambdas = 2)
  test <- evo(X, Y, W, settings)
  
  expect_equal(sum(test$selected), ncol(X))
  
  settings <- get_evo_settings()
  test <- evo(X, Y, W, settings)
  
  expect_equal(sum(test$selected[sol]), pCausal)
  expect_equal(test$c[sol], rep(96.15385, pCausal), tolerance = 1e-5)
  
})