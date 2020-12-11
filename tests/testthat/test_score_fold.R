testthat::skip('TODO')
folds <- lapply(scores, run_scones, eta, lambda, net_matrix)
folds <- do.call(rbind, folds)

criterion <- NULL
K <- NULL
gwas <- NULL
covars <- NULL
net <- NULL

score_fold <- function(folds, criterion, K, gwas, covars, net)

martini:::score_fold(folds, criterion, K, minigwas, covars, net)


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
  scores <- c(95.4, 95.4, 94.1, 96.1, 94.1, 93.3, 0.0, 0.3, 
              0.3, 0.3, 0.3, 0.0, 0.0, 0.0, 96.1, 96.1, 95.1, 
              95.4, 95.4, 96.1, 0.2, 0.3, 0.0, 0.3, 0.0)
  
  skip_on_os("windows")
  expect_equal(cones$selected, cones$c > 0)
  expect_equal(cones$c, scores, tolerance = .1)
  
})