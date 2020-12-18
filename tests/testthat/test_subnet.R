test_that("we recover all the SNPs we want", {
  
  subgi <- martini:::subnet(test_gi, "gene", "A")
  snpA <- grep("A", test_gwas$map$snp.names, value = TRUE)
  
  expect_true(all(names(igraph::V(subgi)) %in% snpA))
  expect_true(all(snpA %in% names(igraph::V(subgi))))
  
})

test_that("we don't recover the SNPs we don't want", {
  
  subgi <- martini:::subnet(test_gi, "gene", "A", affirmative = FALSE)
  snpNotA <- grep("A", test_gwas$map$snp.names, value = TRUE, invert = TRUE)
  
  expect_true(all(names(igraph::V(subgi)) %in% snpNotA))
  expect_true(all(snpNotA %in% names(igraph::V(subgi))))
  
})
