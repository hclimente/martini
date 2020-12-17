testthat::skip_on_bioc()
skip_if(skip_long)

# SNPs from https://www.ebi.ac.uk/gwas/search?query=breast%20cancer
# mapped to fgfr2, tox3, fgfr2, none and none, respectively
brca <- list()
brca$map <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2
                       10 rs2981579 0 121577821 A G
                       16 rs3803662 0 52552429 A G
                       10 rs2981582 0 121592803 A G
                       11 rs614367 0 69513996 A G
                       2 rs13387042 0 217041109 A G
                       ", header = TRUE, stringsAsFactors = FALSE)
brca$map$gpos <- as.numeric(brca$map$gpos)
brca_mapped <- martini:::snp2ensembl(brca)

# 7 snps from https://easygwas.ethz.ch/gwas/results/manhattan/view/4d00706f-ad0f-4f57-9f4e-ac3099b15b94/
# the last one was not assigned to any gene
athal <- list()
athal$map <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2
                       4 Chr4_8297892 0 8297892 A G
                       4 Chr4_8297535 0 8297535 A G
                       4 Chr4_8301059 0 8301059 A G
                       4 Chr4_8300836 0 8300836 A G
                       4 Chr4_8254521 0 8254521 A G
                       4 Chr4_8274507 0 8274507 A G
                       5 Chr5_6485290 0 6485290 A G
                       ", header = TRUE, stringsAsFactors = FALSE)

athal$map$gpos <- as.numeric(athal$map$gpos)
athal_mapped <- martini:::snp2ensembl(athal, organism=3702)

test_that("output is as expected", {
  # dimensions
  expect_equal(ncol(brca_mapped), 2)

  # column order
  expect_equal(length(intersect(brca$map$snp.names, brca_mapped[,1])), 2)
  expect_equal(length(unique(brca_mapped[,2])), 1)
})

test_that("we map snps to their known genes", {

  # correct mapping in humans
  expect_equal(brca_mapped$gene[brca_mapped$snp == "rs2981579"], "FGFR2")
  expect_equal(brca_mapped$gene[brca_mapped$snp == "rs2981582"], "FGFR2")
  expect_equal(length(brca_mapped$gene[brca_mapped$snp == "rs13387042"]), 0)

  # correct mapping in arabidopsis
  expect_equal(athal_mapped$gene[athal_mapped$snp == "Chr4_8297892"], "ACD6")
  expect_equal(athal_mapped$gene[athal_mapped$snp == "Chr4_8297535"], "ACD6")
  expect_equal(athal_mapped$gene[athal_mapped$snp == "Chr4_8301059"], "BHLH104")
  expect_equal(athal_mapped$gene[athal_mapped$snp == "Chr4_8300836"], "BHLH104")
  expect_equal(athal_mapped$gene[athal_mapped$snp == "Chr4_8254521"], "AT4G14342")
  expect_equal(athal_mapped$gene[athal_mapped$snp == "Chr4_8274507"], "AT4G14368")
  expect_equal(length(athal_mapped$gene[athal_mapped$snp == "Chr5_6485290"]), 0)

  # check flanks
  # off by one base to FGFR2 10: 121,478,332-121,598,458 
  red <- list()
  red$map <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2
                       10 outbound 0 121478331 A G
                       ", header = TRUE, stringsAsFactors = FALSE)
  red$map$gpos <- as.numeric(red$map$gpos)

  expect_equal(nrow(martini:::snp2ensembl(red, flank = 1)), 1)
  expect_equal(nrow(martini:::snp2ensembl(red)), NULL)

})
