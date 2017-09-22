library(martini)

snps <- data.frame(snp1 = c("rs1165182","rs9366633","rs1324088","rs9366633","rs1165182","rs9366633"),
                   snp2 = c("rs1408273","rs6905614","rs1324087","rs72843569","rs4712976","rs9393672"))
net <- igraph::graph_from_data_frame(snps, directed = FALSE)

ldnet <- ldweight_edges(net)

test_that("snps in perfect ld are disconnected", {
  expect_false(igraph::are_adjacent(ldnet, "rs1165182", "rs1408273"))
  expect_false(igraph::are_adjacent(ldnet, "rs1324087", "rs1324088"))
  expect_equal(length(igraph::E(net)), length(igraph::E(ldnet)) + 2)
})

test_that("snps not in perfect ld connected", {
  expect_true(all(igraph::E(ldnet) %in% igraph::E(net)))
})
