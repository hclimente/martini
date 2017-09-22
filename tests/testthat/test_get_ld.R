library(martini)

snps <- data.frame(snp1 = c("rs1165182","rs9366633","rs1324088","rs9366633","rs1165182","rs9366633"),
                   snp2 = c("rs1408273","rs6905614","rs1324087","rs72843569","rs4712976","rs9393672"))
net <- igraph::graph_from_data_frame(snps, directed = FALSE)

ld <- get_ld(net)

test_that("output has the right shape", {
  expect_equal(ncol(ld), 2)
  expect_true("key" %in% colnames(ld))
  expect_true("r2" %in% colnames(ld))
})

test_that("we get complete LD-info about the network", {
  expect_equal(length(E(net)), nrow(ld))
  expect_equal(nrow(snps), nrow(ld))
})

test_that("we retrieve known LD values", {
  expect_equal(ld$r2[ld$key == "rs1165182_rs1408273"], 1)
  expect_equal(ld$r2[ld$key == "rs1324087_rs1324088"], 1)
  expect_equal(ld$r2[ld$key == "rs1165182_rs4712976"], 0.548716)
  expect_equal(ld$r2[ld$key == "rs72843569_rs9366633"], 0.069162)
})