library(martini)
source("minimum_data.R")

test_that("we connect the consecutive SNPs", {
  expect_true(are_adjacent(gs, "rs1", "rs2"))
  expect_true(are_adjacent(gs, "rs2", "rs3"))
  expect_false(are_adjacent(gs, "rs3", "rs4"))
})

test_that("we add genomic information to the vertices", {
  expect_equal(get.vertex.attribute(gs, "chr", match(c("rs1", "rs2", "rs3"), V(gs)$name) ), rep(1, 3) )
  expect_equal(get.vertex.attribute(gs, "chr", match(c("rs4", "rs5", "rs6"), V(gs)$name) ), rep(2, 3) )
  expect_equal(get.vertex.attribute(gs, "pos", match(c("rs1", "rs2", "rs3"), V(gs)$name) ), c(10, 20, 30) )
  expect_equal(get.vertex.attribute(gs, "pos", match(c("rs4", "rs5", "rs6"), V(gs)$name) ), c(15, 25, 35) )
})