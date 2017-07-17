library(martini)

gwas <- list()
gwas$map <- data.frame(chr = c(1, 1, 1, 2, 2, 2),
                       snp.names = paste0("rs", 1:6),
                       cm = rep(0, 6),
                       gpos = c(10, 20, 30, 15, 25, 35),
                       allele.1 = rep("A", 6),
                       allele.2 = rep("T", 6))
snp2gene <- data.frame(snp = c("rs1", "rs2", "rs5", "rs6"),
                       gene = c("A", "A", "B", "B"), 
                       stringsAsFactors = FALSE)
ppi <- data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)

net <- get_GI_network(gwas, snp2gene, ppi)
causal <- simulate_causal_snps(net, 3)

test_that("we get causal SNPs from different genes", {
  expect_true((causal$gene %>% unique %>% length) > 1)
  expect_true(are_adjacent(net, causal[1], causal[2]))
  expect_true(are_adjacent(net, causal[2], causal[3]))
})
