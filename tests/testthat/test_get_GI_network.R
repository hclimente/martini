library(martini)
source("minimum_data.R")

test_that("we interconnect the right genes", {
  expect_true(are_adjacent(gi, "rs1", "rs5"))
  expect_false(
    are_adjacent(
      get_GI_network(gwas, 
                     data.frame(snp = c("rs1", "rs2", "rs5", "rs6"),
                                gene = c("A", "A", "C", "B"), 
                                stringsAsFactors = FALSE),
                     data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)),
      "rs1", "rs5"))
})

test_that("we add genomic information to the vertices", {
  expect_equal(get.vertex.attribute(gi, "chr", match(c("rs1", "rs2", "rs3"), V(gi)$name) ), rep(1, 3) )
  expect_equal(get.vertex.attribute(gi, "chr", match(c("rs4", "rs5", "rs6"), V(gi)$name) ), rep(2, 3) )
  expect_equal(get.vertex.attribute(gi, "pos", match(c("rs1", "rs2", "rs3"), V(gi)$name) ), c(10, 20, 30) )
  expect_equal(get.vertex.attribute(gi, "pos", match(c("rs4", "rs5", "rs6"), V(gi)$name) ), c(15, 25, 35) )
})

test_that("we add gene information to the vertices", {
  expect_equal(get.vertex.attribute(gi, "gene", match(c("rs1", "rs2"), V(gi)$name) ), rep("A", 2) )
  expect_equal(get.vertex.attribute(gi, "gene", match(c("rs5", "rs6"), V(gi)$name) ), rep("B", 2) )
  expect_output(get.vertex.attribute(gi, "gene", match("rs3", V(gi)$name) ), NA )
  expect_output(get.vertex.attribute(gi, "gene", match("rs4", V(gi)$name) ), NA )
})

test_that("we are simplifying the network", { 
  
  s2g <- data.frame(snp = c("rs1", "rs2", "rs3", "rs4"),
                    gene = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
  p <- data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)
  x <- as_adj(get_GI_network(gwas, s2g, p))
  
  expect_equal(sum(x != 0 & x != 1), 0)
  })

test_that("warns if ppi is insufficient to create a GI network", {
  expect_warning(get_GI_network(gwas, snp2gene, data.frame(gene1 = "A", gene2 = "A", stringsAsFactors = FALSE)), 
                 "no matches between genes in snp2gene and PPI. No information about PPI will be added.")
  expect_warning(get_GI_network(gwas, snp2gene, data.frame(gene1 = c("A", "B"), gene2 = c("A", "B"), stringsAsFactors = FALSE)), 
                 "no matches between genes in snp2gene and PPI. No information about PPI will be added.")
})