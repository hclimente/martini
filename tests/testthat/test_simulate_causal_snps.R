library(martini)
source("minimum_data.R")

# GI network
causal2 <- simulate_causal_snps(gi, 2)
causal3 <- simulate_causal_snps(gi, 3)

test_that("we get causal SNPs from two different genes", {
  expect_equal(length(causal2$gene %>% unique), 2)
  expect_equal(length(causal3$gene %>% unique), 3)
})

test_that("genes with less than 2 single-gene SNPs are discarded", {
  expect_false("D" %in% causal2$gene)
  expect_false("D" %in% causal3$gene)
})

test_that("SNPs are interconnected", {
  expect_equal((igraph::components(igraph::induced_subgraph(gi, names(causal2))) %>% 
                .$membership %>% 
                unique), 1 )
  expect_equal((igraph::components(igraph::induced_subgraph(gi, names(causal3))) %>% 
                  .$membership %>% 
                  unique), 1 )
})

half <- simulate_causal_snps(gi, 2, 0.5)
test_that("we can modulate the proportion of SNPs", {
  expect_true(length(causal2) == (length(half) * 2))
})