library(martini)

ld <- snpStats::ld(minigwas$genotypes, depth = 2, stats = "R.squared")
pruned <- ld_prune(minigwas, ld)

allLd <- ld
allLd[is.na(ld)] <- 1
allPruned <- ld_prune(minigwas, allLd)

le <- ld
le[ld == 1] <- NA
zeroPruned <- ld_prune(minigwas, le)

test_that("output is as expected", {
  
  # no changes to individual information
  expect_equal(pruned$fam, minigwas$fam)
  
  # appropriate changes to map
  expect_true(nrow(pruned$map) < nrow(minigwas$map))
  expect_equal(ncol(pruned$map), ncol(minigwas$map))
  expect_true(is.data.frame(pruned$map))
  
  # appropriate changes to genotypes
  expect_equal(nrow(pruned$genotypes), nrow(minigwas$genotypes))
  expect_true(ncol(pruned$genotypes) < ncol(minigwas$genotypes))
  expect_equal(class(pruned$genotypes)[1], "SnpMatrix")
  expect_equal(colnames(pruned$genotypes), pruned$map$snp.names)
  
  # prune information
  expect_true("pruning" %in% names(pruned))
  
})

test_that("it prunes", {
  
  # we dont remove SNPs where we have NAs
  expect_equal(sum(grepl("[BD-]", pruned$map$snp.names)), 
               sum(grepl("[BD-]", minigwas$map$snp.names)))
  expect_equal(length(unique(pruned$pruning$block[grepl("[BD-]", pruned$pruning$snp)])), 13)
  
  # we remove SNPs where correlations are perfect i.e. same haploblock
  expect_equal(sum(grepl("[AC]", pruned$map$snp.names)) + 10, 
               sum(grepl("[AC]", minigwas$map$snp.names)))
  expect_equal(length(unique(pruned$pruning$block[grepl("[AC]", pruned$pruning$snp)])), 2)
  
  # only remove SNPs that we are supposed to
  expect_equal(nrow(pruned$map) + 10, nrow(minigwas$map))
  
  # extreme cases
  expect_equal(nrow(allPruned$map), 1)
  expect_equal(zeroPruned$map, minigwas$map)
  
})

test_that("cutoff works", {
  
  test <- ld_prune(minigwas, ld, cutoff = 2)
  
  # we dont remove any SNPs where we have NAs
  expect_equal(test$map, minigwas$map)
  expect_equal(test$genotypes, minigwas$genotypes)
  
})

test_that("SNP order is checked", {
  
  shuffledSNPs <- sample(rownames(ld), 25)
  shuffledLD <- ld[shuffledSNPs, shuffledSNPs]
  
  expect_error(ld_prune(minigwas, shuffledLD), "gwas$map and ld SNP order differ.", fixed = TRUE)
  
})