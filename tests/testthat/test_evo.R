library(martini)

miniX <- as(minigwas$genotypes, "numeric")
miniY <- minigwas$fam$affected
gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
miniW <- igraph::as_adj(gi)
miniW <- miniW[colnames(miniX), colnames(miniX)]

test_that("output is as expected", {
  
  settings <- get_evo_settings(etas = c(1,2,3), lambdas = c(2.5,3.5))
  test <- evo(miniX, miniY, miniW, settings)
  
  expect_equal(length(test), 7)
  expect_equal(class(test), "list")
  expect_true(test$eta %in% settings[['etas']])
  expect_true(test$lambda %in% settings[['lambdas']])
  expect_equal(dim(test$grid), c(3,2))
  expect_equal(class(test$selected), "numeric")
  expect_equal(class(test$c), "numeric")
  
})

test_that("we recover causal SNPs", {
  
  settings <- get_evo_settings(etas = 1, lambdas = 2)
  test <- evo(miniX, miniY, miniW, settings)
  
  expect_equal(sum(test$selected), ncol(miniX))
  
  settings <- get_evo_settings()
  test <- evo(miniX, miniY, miniW, settings)
  
  skip_on_os("windows")
  expect_equal(test$selected, as.numeric(test$c > 0))
  expect_equal(test$c[test$selected], rep(95.5, sum(test$selected)), tolerance = .1)
  
})