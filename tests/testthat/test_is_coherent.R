library(martini)

suffledMap <- minigwas$map[sample(1:nrow(minigwas$map)), ]
shuffledGenotypes <- minigwas$genotypes[, suffledMap$snp.names]

test_that("we get errors when we should", {
  
  expect_true(is_coherent(minigwas))
  expect_error(is_coherent(list(map = suffledMap,
                                genotypes = minigwas$genotypes,
                                fam = minigwas$fam)),
               "map is not ordered by genomic position.")
  expect_error(is_coherent(list(map = minigwas$map,
                                genotypes = shuffledGenotypes,
                                fam = minigwas$fam)),
               "map and genotype SNP order differ.")
  
})
