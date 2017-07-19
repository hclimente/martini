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