library(martini)

# get_ld_from_net
snps <- data.frame(snp1 = c("rs1165182","rs9366633","rs1324088","rs9366633","rs1165182","rs9366633"),
                   snp2 = c("rs1408273","rs6905614","rs1324087","rs72843569","rs4712976","rs9393672"))
net <- igraph::graph_from_data_frame(snps, directed = FALSE)
ld <- get_ld(net)

test_that("get_ld_from_net output has the right shape", {
  expect_equal(ncol(ld), 2)
  expect_true("key" %in% colnames(ld))
  expect_true("r2" %in% colnames(ld))
})

test_that("get_ld_from_net returns all LD-info available in Ensembl", {
  # rs9366633--rs72843569 LD unavailable for CEU
  expect_equal(nrow(ld), length(E(net)) - 1)
})

test_that("get_ld_from_net returns known LD values", {
  expect_equal(ld$r2[ld$key == "rs1165182_rs1408273"], 0.98)
  expect_equal(ld$r2[ld$key == "rs1324087_rs1324088"], 1)
  expect_equal(ld$r2[ld$key == "rs1165182_rs4712976"], 0.296567)
  expect_equal(ld$r2[ld$key == "rs9366633_rs9393672"], 0.319992)
})


# get_ld_from_gwas
gwas <- list()
gwas$map <- data.frame(V1 = c(9,9,9,9,9),
                         V2 = c("rs1333049","rs7857118","rs1333047","rs7862252","rs7848340"),
                       V3 = 0,
                       V4 = c(22125504,22124141,22124505,24125064,24125020), 
                       V5 = "A", V6 = "C")
ld <- get_ld(gwas = gwas)

test_that("get_ld_from_gwas output has the right shape", {
  expect_equal(ncol(ld), 2)
  expect_true("key" %in% colnames(ld))
  expect_true("r2" %in% colnames(ld))
})

test_that("get_ld_from_gwas returns all LD-info available in Ensembl", {
  expect_equal(nrow(ld), 4)
})

test_that("get_ld_from_gwas returns known LD values", {
  expect_equal(ld$r2[ld$key == "rs1333049_rs7857118"], 0.882052)
  expect_equal(ld$r2[ld$key == "rs1333047_rs1333049"], 0.959308)
  # further away than 1Mb
  expect_false("rs1333047_rs7848340" %in% ld$key)
})