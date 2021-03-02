test_that("we connect the consecutive SNPs", {
  expect_true(igraph::are_adjacent(mini_gs, "rs1", "rs2"))
  expect_true(igraph::are_adjacent(mini_gs, "rs2", "rs3"))
  expect_false(igraph::are_adjacent(mini_gs, "rs3", "rs4"))
})

test_that("we add genomic information to the vertices", {
  expect_equal(igraph::get.vertex.attribute(mini_gs, "chr", match(c("rs1", "rs2", "rs3"), igraph::V(mini_gs)$name) ), rep(1, 3) )
  expect_equal(igraph::get.vertex.attribute(mini_gs, "chr", match(c("rs4", "rs5", "rs6"), igraph::V(mini_gs)$name) ), rep(2, 3) )
  expect_equal(igraph::get.vertex.attribute(mini_gs, "pos", match(c("rs1", "rs2", "rs3"), igraph::V(mini_gs)$name) ), c(10, 20, 30) )
  expect_equal(igraph::get.vertex.attribute(mini_gs, "pos", match(c("rs4", "rs5", "rs6"), igraph::V(mini_gs)$name) ), c(15, 25, 35) )
})

test_that("edges have weights", {
  expect_false(is.null(igraph::E(mini_gs)$weight))
})
