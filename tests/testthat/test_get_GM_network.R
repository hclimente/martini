library(martini)

gwas <- list()
gwas$map <- data.frame(chr = c(1, 1, 1, 2, 2, 2),
                       snp.names = paste0("rs", 1:6),
                       cm = rep(0, 6),
                       gpos = c(10, 20, 30, 15, 25, 35),
                       allele.1 = rep("A", 6),
                       allele.2 = rep("T", 6))

test_that("we interconnect snps from a gene", {
  expect_true(are_adjacent(get_GM_network(gwas, data.frame(snp = c("rs1", "rs2"), gene = "1")), "rs1", "rs2"))
  expect_true(are_adjacent(get_GM_network(gwas, data.frame(snp = c("rs1", "rs2", "rs3"), gene = "1")), "rs1", "rs3"))
  expect_false(are_adjacent(get_GM_network(gwas, data.frame(snp = c("rs1", "rs2", "rs6"), gene = c("1","1","2"))), "rs1", "rs6"))
})

test_that("crash if snp2gene is insufficient to create a GM network", {
  expect_error(get_GM_network(gwas, data.frame(snp = "rs1", gene = "A")), "the data frame should contain at least two columns")
})

net <- get_GM_network(gwas, data.frame(snp = c("rs1", "rs2", "rs5", "rs6"), 
                                       gene = c("A", "A", "B", "B"), 
                                       stringsAsFactors = FALSE))

test_that("we add genomic information to the vertices", {
  expect_equal(get.vertex.attribute(net, "chr", match(c("rs1", "rs2", "rs3"), V(net)$name) ), rep(1, 3) )
  expect_equal(get.vertex.attribute(net, "chr", match(c("rs4", "rs5", "rs6"), V(net)$name) ), rep(2, 3) )
  expect_equal(get.vertex.attribute(net, "pos", match(c("rs1", "rs2", "rs3"), V(net)$name) ), c(10, 20, 30) )
  expect_equal(get.vertex.attribute(net, "pos", match(c("rs4", "rs5", "rs6"), V(net)$name) ), c(15, 25, 35) )
})

test_that("we add gene information to the vertices", {
  expect_equal(get.vertex.attribute(net, "gene", match(c("rs1", "rs2"), V(net)$name) ), rep("A", 2) )
  expect_equal(get.vertex.attribute(net, "gene", match(c("rs5", "rs6"), V(net)$name) ), rep("B", 2) )
  expect_output(get.vertex.attribute(net, "gene", match("rs3", V(net)$name) ), NA )
  expect_output(get.vertex.attribute(net, "gene", match("rs4", V(net)$name) ), NA )
})