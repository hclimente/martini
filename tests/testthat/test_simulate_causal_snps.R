library(martini)
source("minimum_data.R")

net <- get_GI_network(gwas, snp2gene, ppi)
causal <- simulate_causal_snps(net, 3)

test_that("we get causal SNPs from different genes", {
  expect_true((causal$gene %>% unique %>% length) > 1)
  expect_true(are_adjacent(net, causal[1], causal[2]))
  expect_true(are_adjacent(net, causal[2], causal[3]))
})
