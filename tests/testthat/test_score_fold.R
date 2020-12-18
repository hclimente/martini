K <- cut(seq(1, nrow(minigwas[['fam']])), breaks = 2, labels = FALSE)

test_that("stability works", {
  
  criterion <- 'stability'
  
  selected <- c(T,T,T,T,T,F,F,F,F,F)
  expect_equal(selected, martini:::score_fold(minigwas, data.frame(), gi, selected, criterion))
  
})

test_that("clustering coefficients works", {
  
  # test edge cases 
  full_graph <- make_full_graph(10)
  test_graph <- full_graph + graph_from_edgelist(cbind(seq(1, 14), seq(2, 15)), directed = FALSE)
  V(test_graph)$name <- minigwas[['map']][['snp.names']]
  
  selected <- c(rep(TRUE, 10), rep(FALSE, 15))
  expect_equal(1, martini:::score_fold(minigwas, data.frame(), test_graph, selected, 'global_clustering'))
  expect_equal(1, martini:::score_fold(minigwas, data.frame(), test_graph, selected, 'local_clustering'))
  
  selected <- c(rep(FALSE, 15), rep(TRUE, 10))
  expect_equal(0, martini:::score_fold(minigwas, data.frame(), test_graph, selected, 'global_clustering'))
  expect_equal(0, martini:::score_fold(minigwas, data.frame(), test_graph, selected, 'local_clustering'))
  
})
