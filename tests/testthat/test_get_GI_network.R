source("minimum_data.R")

test_that("we interconnect the right genes", {
  expect_true(igraph::are_adjacent(gi, "rs1", "rs5"))
  expect_false(
    igraph::are_adjacent(
      get_GI_network(gwas, 
                     snpMapping = data.frame(snp = c("rs1", "rs2", "rs5", "rs6"),
                                             gene = c("A", "A", "C", "B"), 
                                             stringsAsFactors = FALSE),
                     ppi = data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)),
      "rs1", "rs5"))
})

test_that("we add genomic information to the vertices", {
  expect_equal(igraph::get.vertex.attribute(gi, "chr", match(c("rs1", "rs2", "rs3"), igraph::V(gi)$name) ), rep(1, 3) )
  expect_equal(igraph::get.vertex.attribute(gi, "chr", match(c("rs4", "rs5", "rs6"), igraph::V(gi)$name) ), rep(2, 3) )
  expect_equal(igraph::get.vertex.attribute(gi, "pos", match(c("rs1", "rs2", "rs3"), igraph::V(gi)$name) ), c(10, 20, 30) )
  expect_equal(igraph::get.vertex.attribute(gi, "pos", match(c("rs4", "rs5", "rs6"), igraph::V(gi)$name) ), c(15, 25, 35) )
})

test_that("we add gene information to the vertices", {
  expect_equal(igraph::get.vertex.attribute(gi, "gene", match(c("rs1", "rs2"), igraph::V(gi)$name) ), rep("A", 2) )
  expect_equal(igraph::get.vertex.attribute(gi, "gene", match(c("rs5", "rs6"), igraph::V(gi)$name) ), rep("B", 2) )
  expect_output(igraph::get.vertex.attribute(gi, "gene", match("rs3", igraph::V(gi)$name) ), NA )
  expect_output(igraph::get.vertex.attribute(gi, "gene", match("rs4", igraph::V(gi)$name) ), NA )
})

test_that("we are simplifying the network", { 
  
  s2g <- data.frame(snp = c("rs1", "rs2", "rs3", "rs4"),
                    gene = c("A", "A", "B", "B"), stringsAsFactors = FALSE)
  p <- data.frame(gene1 = "A", gene2 = "B", stringsAsFactors = FALSE)
  x <- igraph::as_adj(get_GI_network(gwas, snpMapping = s2g, ppi = p))
  
  expect_equal(sum(x != 0 & x != 1), 0)
  })

test_that("warns if ppi is insufficient to create a GI network", {
  expect_warning(get_GI_network(gwas, 
                                snpMapping = snpMapping, 
                                ppi = data.frame(gene1 = "A", gene2 = "A", stringsAsFactors = FALSE)), 
                 "no matches between genes in snpMapping and PPI. No information about PPI will be added.")
  expect_warning(get_GI_network(gwas, 
                                snpMapping = snpMapping, 
                                ppi = data.frame(gene1 = c("A", "B"), gene2 = c("A", "B"), stringsAsFactors = FALSE)), 
                 "no matches between genes in snpMapping and PPI. No information about PPI will be added.")
})

test_that("edges have weights", {
  expect_false(is.null(igraph::E(gi)$weight))
})