library(martini)
source("minimum_data.R")

test_that("we interconnect snps from a gene", {
  expect_true(igraph::are_adjacent(get_GM_network(gwas, snpMapping = data.frame(snp = c("rs1", "rs2"), gene = "1")), "rs1", "rs2"))
  expect_true(igraph::are_adjacent(get_GM_network(gwas, snpMapping = data.frame(snp = c("rs1", "rs2", "rs3"), gene = "1")), "rs1", "rs3"))
  expect_false(igraph::are_adjacent(get_GM_network(gwas, snpMapping = data.frame(snp = c("rs1", "rs2", "rs6"), gene = c("1","1","2"))), "rs1", "rs6"))
})

test_that("warns if snpMapping is insufficient to create a GM network", {
  expect_warning(get_GM_network(gwas, snpMapping = data.frame(snp = "rs1", gene = "A")), 
                 "insufficient information to add gene information")
  expect_warning(get_GM_network(gwas, snpMapping = data.frame(snp = c("rs1", "rs2"), gene = c("A", "B"))), 
                 "insufficient information to add gene information")
})

test_that("we add genomic information to the vertices", {
  expect_equal(igraph::get.vertex.attribute(gm, "chr", match(c("rs1", "rs2", "rs3"), igraph::V(gm)$name) ), rep(1, 3) )
  expect_equal(igraph::get.vertex.attribute(gm, "chr", match(c("rs4", "rs5", "rs6"), igraph::V(gm)$name) ), rep(2, 3) )
  expect_equal(igraph::get.vertex.attribute(gm, "pos", match(c("rs1", "rs2", "rs3"), igraph::V(gm)$name) ), c(10, 20, 30) )
  expect_equal(igraph::get.vertex.attribute(gm, "pos", match(c("rs4", "rs5", "rs6"), igraph::V(gm)$name) ), c(15, 25, 35) )
})

test_that("we add gene information to the vertices", {
  expect_equal(igraph::get.vertex.attribute(gm, "gene", match(c("rs1", "rs2"), igraph::V(gm)$name) ), rep("A", 2) )
  expect_equal(igraph::get.vertex.attribute(gm, "gene", match(c("rs5", "rs6"), igraph::V(gm)$name) ), rep("B", 2) )
  expect_output(igraph::get.vertex.attribute(gm, "gene", match("rs3", igraph::V(gm)$name) ), NA )
  expect_output(igraph::get.vertex.attribute(gm, "gene", match("rs4", igraph::V(gm)$name) ), NA )
})

test_that("we are simplifying the network", { 
  
  s2g <- data.frame(snp = c("rs1", "rs2", "rs3", "rs4"),
                    gene = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
  x <- igraph::as_adj(get_GM_network(gwas, snpMapping = s2g))
  
  expect_equal(sum(x != 0 & x != 1), 0)
})