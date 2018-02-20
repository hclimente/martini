gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)

test_that("we recover all the SNPs we want", {
  
  subgi <- martini:::subvert(gi, "gene", "A")
  snpA <- grep("A", minigwas$map$snp.names, value = TRUE)
  
  expect_true(all(names(subgi) %in% snpA))
  expect_true(all(snpA %in% names(subgi)))
  
})

test_that("we don't recover the SNPs we don't want", {
  
  subgi <- martini:::subvert(gi, "gene", "A", affirmative = FALSE)
  snpNotA <- grep("A", minigwas$map$snp.names, value = TRUE, invert = TRUE)
  
  expect_true(all(names(subgi) %in% snpNotA))
  expect_true(all(snpNotA %in% names(subgi)))
  
})
