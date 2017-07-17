library(martini)

gwas <- list()
gwas$map <- data.frame(chr = c(1, 1, 1, 2, 2, 2),
                       snp.names = paste0("rs", 1:6),
                       cm = rep(0, 6),
                       gpos = c(10, 20, 30, 15, 25, 35),
                       allele.1 = rep("A", 6),
                       allele.2 = rep("T", 6))

net <- get_GS_network(gwas)

test_that("we connect the consecutive SNPs", {
  expect_true(are_adjacent(net, "rs1", "rs2"))
  expect_true(are_adjacent(net, "rs2", "rs3"))
  expect_false(are_adjacent(net, "rs3", "rs4"))
})

test_that("we add genomic information to the vertices", {
  expect_equal(get.vertex.attribute(net, "chr", match(c("rs1", "rs2", "rs3"), V(net)$name) ), rep(1, 3) )
  expect_equal(get.vertex.attribute(net, "chr", match(c("rs4", "rs5", "rs6"), V(net)$name) ), rep(2, 3) )
  expect_equal(get.vertex.attribute(net, "pos", match(c("rs1", "rs2", "rs3"), V(net)$name) ), c(10, 20, 30) )
  expect_equal(get.vertex.attribute(net, "pos", match(c("rs4", "rs5", "rs6"), V(net)$name) ), c(15, 25, 35) )
})