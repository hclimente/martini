library(martini)

test_that("output is as expected", {
  
  c <- rep(1, 3)
  W <- matrix(c(0,1,0,1,0,1,0,1,0,0,0,1), 3, 3)
  W <- as(W, "sparseMatrix")
  
  selected <- mincut_c(c, 1, 2, W)
  
  expect_equal(length(selected), 3)
  expect_equal(class(selected), "logical")

})

test_that("we recover causal SNPs", {
  
  pCausal <- 20
  pNonCausal <- 100
  
  c <- c(rep(10, pCausal), rep(0,pNonCausal))
  
  W <- diag(1, pCausal + pNonCausal, pCausal + pNonCausal)
  W <- cbind(rbind(0, W[-1,-1]),0) + rbind(cbind(0,W[-1,-1]),0)
  W[1:pCausal, 1:pCausal] <- 1
  diag(W) <- 0
  W <- as(W, "sparseMatrix")
  
  selected <- mincut_c(c, 1, 2, W)
  
  expect_equal(sum(selected), pCausal)
  expect_true(all(as.logical(selected[1:pCausal])))
  
})
