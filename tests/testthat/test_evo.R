library(martini)

test_that("output is as expected", {
  
  X <- matrix(c(0,0,0,0,1,1,1,1,2,2,2,2), 4, 3)
  Y <- c(1,1,2,2)
  W <- matrix(c(0,1,0,1,0,1,0,1,0,0,0,1), 3, 3)
  W <- as(W, "sparseMatrix")
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
  
  N <- 100
  pCausal <- 20
  pNonCausal <- 100
  causal <- c(rep(2, N/2), rep(0, N/2))
  rest <- rep(0,N)
  
  X <- cbind(matrix(rep(causal, pCausal), N, pCausal),
             matrix(rep(rest, pNonCausal), N, pNonCausal))
  Y <- c(rep(2, N/2), rep(1, N/2))
  
  W <- diag(1, pCausal + pNonCausal, pCausal + pNonCausal)
  W <- cbind(rbind(0, W[-1,-1]),0) + rbind(cbind(0,W[-1,-1]),0)
  W[1:pCausal, 1:pCausal] <- 1
  diag(W) <- 0
  W <- as(W, "sparseMatrix")
  settings <- get_evo_settings(etas = 1, lambdas = 2)
  
  test <- evo(X, Y, W, settings)
  
  expect_equal(sum(test$selected), pCausal)
  expect_true(all(as.logical(test$selected[1:pCausal])))
  expect_equal(test$c[1:pCausal], rep(96.15385, pCausal), tolerance = 1e-5)
  
  settings <- get_evo_settings()
  
  test <- evo(X, Y, W, settings)
  
  expect_equal(sum(test$selected), pCausal)
  expect_true(all(as.logical(test$selected[1:pCausal])))
  expect_equal(test$c[1:pCausal], rep(96.15385, pCausal), tolerance = 1e-5)
  
})