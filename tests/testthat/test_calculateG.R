library(martini)

test_that("effect size modulates the effect", {
  
  X <- matrix(c(0,1,2), 3, 1)
  G1 <- calculateG(1, X, "additive")
  G2 <- calculateG(2, X, "additive")
  G3 <- calculateG(3, X, "additive")
  
  expect_true(all(G1 * 2 == G2))
  expect_true(all(G1 * 3 == G3))
  
})

test_that("SNP frequency modulates the effect", {
  
  X <- matrix(c(0,1,2), 3, 1)
  G <- calculateG(1, X, "additive")
  
  # check that samples with the same number of minor alleles have the same effect
  expect_equal(G[1], -G[3])
  expect_equal(G[2], 0)

  sorted <- sort(rowSums(X), decreasing = TRUE, index.return=TRUE)$ix

  # check that more SNPs result in stronger effects
  expect_false(is.unsorted(G[sorted]))
  
})

test_that("additive effects are the sum of the effects of individual SNPs", {
  
  X1 <- matrix(c(0,1,2), 3, 1)
  X2 <- matrix(c(1,1,2), 3, 1)
  X3 <- matrix(c(0,2,2), 3, 1)
  X <- cbind(X1, X2, X3)
  
  u <- c(1, 3.5, 2.1)
  
  G1 <- calculateG(u[1], X1, "additive")
  G2 <- calculateG(u[2], X2, "additive")
  G3 <- calculateG(u[3], X3, "additive")
  G <- calculateG(u, X, "additive")
  
  expect_equal(G, G1 + G2 + G3)
  
})

test_that("errors when it should", { 
  
  X <- matrix(c(0,1,2), 3, 1)
  
  expect_error(calculateG(1, X, "xyz"), "Genetic model xyz not recognised.", fixed=TRUE)
  
})