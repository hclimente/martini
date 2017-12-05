library(martini)
load("examplegwas.rda")

gwas <- examplegwas$gwas
ld <- data.frame(key = c("rs1101-rs1101", "rs1101-rs1103", "rs1102-rs1103"),
                 r2 = c(0, 1, 0.5), stringsAsFactors = F)

snps <- data.frame(snp1 = c("rs1101","rs1101","rs1102"),
                   snp2 = c("rs1102","rs1103","rs1103"))
net <- igraph::graph_from_data_frame(snps, directed = FALSE)
net <- igraph::set_edge_attr(net, "weight", value = 1)

ldnet <- ldweight_edges(net, ld, method = "subtraction")

test_that("snps in perfect ld are disconnected", {
  expect_false(igraph::are_adjacent(ldnet, "rs1101", "rs1103"))
  expect_equal(length(igraph::E(ldnet)), 2)
})

test_that("snps not in perfect ld connected", {
  expect_true(all(igraph::E(ldnet) %in% igraph::E(net)))
})
