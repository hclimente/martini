genotypes <- minigwas[['genotypes']]
phenotypes <- minigwas[['fam']][['affected']]
covars <- data.frame(sample = sort(minigwas[['fam']][['member']], decreasing = T), 
                     confounding = c(rep(1, nrow(genotypes)/2), rep(10, nrow(genotypes)/2)),
                     unrelated = c(rep(1, nrow(genotypes))))

test_that("chi2 are correctly calculated", {
  
  true_vals <- rep(0, nrow(minigwas[['map']]))
  true_vals[grepl("[AC]", minigwas[['map']][['snp.names']])] <- 99
  names(true_vals) <- minigwas[['map']][['snp.names']]
  c <- single_snp_association(genotypes, phenotypes, covars, 'chi2')
  
  expect_equal(c, true_vals)
  
  # covariates have no effect
  c <- single_snp_association(genotypes, phenotypes, data.frame(), 'chi2')
  expect_equal(c, true_vals)
  
})

test_that("glm's chi2 are correctly calculated", {
  
  # without covariates
  true_vals <- rep(0, nrow(minigwas[['map']]))
  true_vals[grepl("[AC]", minigwas[['map']][['snp.names']])] <- 100
  names(true_vals) <- minigwas[['map']][['snp.names']]
  
  c <- single_snp_association(genotypes, phenotypes, data.frame(), 'glm')
  
  expect_equal(c, true_vals)
  
  # with explanatory covariates
  true_vals <- rep(0, nrow(minigwas[['map']]))
  names(true_vals) <- minigwas[['map']][['snp.names']]
  c <- single_snp_association(genotypes, phenotypes, covars, 'glm')
  
  expect_equal(c, true_vals)
  
})