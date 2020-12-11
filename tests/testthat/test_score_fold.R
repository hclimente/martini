source("big_network.R")

K <- cut(seq(1, nrow(minigwas[['fam']])), breaks = 2, labels = FALSE)
covars <- data.frame()

test_that("consistency works", {
  
  criterion <- 'consistency'
  
  folds <- rbind(c(rep(TRUE, 5), rep(FALSE, 5)), c(rep(TRUE, 5), rep(FALSE, 5)))
  expect_equal(1, martini:::score_fold(folds, criterion, K, minigwas, covars, gi))
  
  folds <- rbind(c(rep(TRUE, 5), rep(FALSE, 5)), c(rep(FALSE, 5), rep(TRUE, 5)))
  expect_equal(0, martini:::score_fold(folds, criterion, K, minigwas, covars, gi))
  
  folds <- rbind(c(F,T,F,T,F,T,F,T,F,T), c(T,T,T,T,T,F,F,F,F,F))
  expect_equal(.25, martini:::score_fold(folds, criterion, K, minigwas, covars, gi))
  
})

test_that("clustering coefficients works", {
  
  # test edge cases 
  test_graph <- full_graph + graph_from_edgelist(cbind(seq(1, 14), seq(2, 15)), directed = FALSE)
  V(test_graph)$name <- minigwas[['map']][['snp.names']]
  
  folds <- rbind(c(rep(TRUE, 10), rep(FALSE, 15)), c(rep(TRUE, 10), rep(FALSE, 15)))
  expect_equal(2, martini:::score_fold(folds, 'global_clustering', K, minigwas, covars, test_graph, 1))
  expect_equal(2, martini:::score_fold(folds, 'local_clustering', K, minigwas, covars, test_graph, 1))
  
  folds <- rbind(c(rep(FALSE, 15), rep(TRUE, 10)), c(rep(FALSE, 15), rep(TRUE, 10)))
  expect_equal(0, martini:::score_fold(folds, 'global_clustering', K, minigwas, covars, test_graph, 1))
  expect_equal(0, martini:::score_fold(folds, 'local_clustering', K, minigwas, covars, test_graph, 1))
  
  # test empirical cases
  set.seed(0)
  selected <- unlist(scones.cv(minigwas, gi)['selected'])
  selected2 <- sample(selected, length(selected))
  
  folds <- rbind(selected, selected2)
  expect_equal(1.56, martini:::score_fold(folds, 'global_clustering', K, minigwas, covars, gi), tolerance = 0.01)
  expect_equal(1.68, martini:::score_fold(folds, 'local_clustering', K, minigwas, covars, gi), tolerance = 0.01)
  
})
