library(martini)
source("minimum_data.R")

test_that("we interconnect snps from a gene", {
  expect_true(are_adjacent(get_GM_network(gwas, data.frame(snp = c("rs1", "rs2"), gene = "1")), "rs1", "rs2"))
  expect_true(are_adjacent(get_GM_network(gwas, data.frame(snp = c("rs1", "rs2", "rs3"), gene = "1")), "rs1", "rs3"))
  expect_false(are_adjacent(get_GM_network(gwas, data.frame(snp = c("rs1", "rs2", "rs6"), gene = c("1","1","2"))), "rs1", "rs6"))
})

test_that("crash if snp2gene is insufficient to create a GM network", {
  expect_error(get_GM_network(gwas, data.frame(snp = "rs1", gene = "A")), "the data frame should contain at least two columns")
})

test_that("we add genomic information to the vertices", {
  expect_equal(get.vertex.attribute(gm, "chr", match(c("rs1", "rs2", "rs3"), V(gm)$name) ), rep(1, 3) )
  expect_equal(get.vertex.attribute(gm, "chr", match(c("rs4", "rs5", "rs6"), V(gm)$name) ), rep(2, 3) )
  expect_equal(get.vertex.attribute(gm, "pos", match(c("rs1", "rs2", "rs3"), V(gm)$name) ), c(10, 20, 30) )
  expect_equal(get.vertex.attribute(gm, "pos", match(c("rs4", "rs5", "rs6"), V(gm)$name) ), c(15, 25, 35) )
})

test_that("we add gene information to the vertices", {
  expect_equal(get.vertex.attribute(gm, "gene", match(c("rs1", "rs2"), V(gm)$name) ), rep("A", 2) )
  expect_equal(get.vertex.attribute(gm, "gene", match(c("rs5", "rs6"), V(gm)$name) ), rep("B", 2) )
  expect_output(get.vertex.attribute(gm, "gene", match("rs3", V(gm)$name) ), NA )
  expect_output(get.vertex.attribute(gm, "gene", match("rs4", V(gm)$name) ), NA )
})