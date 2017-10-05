library(martini)

gwas <- examplegwas$gwas
ld <- get_ld(gwas)

test_that("get_ld_from_gwas output has the right shape", {
  expect_equal(ncol(ld), 2)
  expect_true("key" %in% colnames(ld))
  expect_true("r2" %in% colnames(ld))
  expect_equal(nrow(ld), nrow(unique(ld)))
})

X <- as(gwas$genotypes, "numeric")
rs275 <- X[gwas$fam$affected == 1, gwas$map$snp.names=="rs275"]
rs276 <- X[gwas$fam$affected == 1, gwas$map$snp.names=="rs276"]
rs1109 <- X[gwas$fam$affected == 1, gwas$map$snp.names=="rs1109"]
rs1110 <- X[gwas$fam$affected == 1, gwas$map$snp.names=="rs1110"]
rs1717 <- X[gwas$fam$affected == 1, gwas$map$snp.names=="rs1717"]
rs1718 <- X[gwas$fam$affected == 1, gwas$map$snp.names=="rs1718"]

test_that("get_ld_from_gwas returns known LD values", {
  expect_equal(ld$r2[ld$key == "rs275-rs276"], cor(rs275, rs276)^2)
  expect_equal(ld$r2[ld$key == "rs1109-rs1110"], cor(rs1109, rs1110)^2)
  expect_equal(ld$r2[ld$key == "rs1717-rs1718"], cor(rs1717, rs1718)^2)
  # further away than 1Mb
  expect_false("rs1101-rs988" %in% ld$key)
})