K <- cut(seq(1, nrow(test_gwas[['fam']])), breaks = 2, labels = FALSE)

test_that("stability works", {
  
  criterion <- 'stability'
  
  selected <- c(T,T,T,T,T,F,F,F,F,F)
  expect_equal(selected, martini:::score_fold(test_gwas, data.frame(), gi, selected, criterion))
  
})

test_that("clustering coefficients works", {
  
  # test edge cases 
  full_graph <- make_full_graph(10)
  test_graph <- full_graph + graph_from_edgelist(cbind(seq(1, 14), seq(2, 15)), directed = FALSE)
  V(test_graph)$name <- test_gwas[['map']][['snp.names']]
  
  selected <- c(rep(TRUE, 10), rep(FALSE, 15))
  expect_equal(1, martini:::score_fold(test_gwas, data.frame(), test_graph, selected, 'global_clustering'))
  expect_equal(1, martini:::score_fold(test_gwas, data.frame(), test_graph, selected, 'local_clustering'))
  
  selected <- c(rep(FALSE, 15), rep(TRUE, 10))
  expect_equal(0, martini:::score_fold(test_gwas, data.frame(), test_graph, selected, 'global_clustering'))
  expect_equal(0, martini:::score_fold(test_gwas, data.frame(), test_graph, selected, 'local_clustering'))
  
})

test_that('log-penalized functions work', {
  
  selected <- c(T,T,T,T,T,F,F,F,F,F)
  phenotypes <- test_gwas[['fam']][['affected']]
  genotypes <- as(test_gwas[['genotypes']], 'numeric')
  genotypes <- as.data.frame(genotypes[, selected])
  genotypes <- cbind(genotypes, covars[,c('confounding','unrelated')])
  
  k <- ncol(genotypes)
  n <- nrow(genotypes)
  
  set.seed(0)
  model <- glm(phenotypes ~ ., data = genotypes)
  expect_equal(gsub('`', '', names(coefficients(model))), 
               c("(Intercept)", colnames(genotypes)))
  
  set.seed(0)
  expect_equal(-BIC(model), martini:::score_fold(test_gwas, covars, gi, selected, 'bic'), tolerance = 4)
  expect_equal(-AIC(model), martini:::score_fold(test_gwas, covars, gi, selected, 'aic'), tolerance = 4)
  expect_equal(-(AIC(model) + (2*k^2 + 2*k)/(n-k-1)), 
               martini:::score_fold(test_gwas, covars, gi, selected, 'aicc'), tolerance = 4)
  
})
