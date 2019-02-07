genotypes <- minigwas[['genotypes']]
phenotypes <- minigwas[['fam']][['affected']]
covars <- data.frame(sample = sort(minigwas[['fam']][['member']], decreasing = T), 
                     confounding = c(rep(1, nrow(genotypes)/2), rep(10, nrow(genotypes)/2)),
                     unrelated = c(rep(1, nrow(genotypes))))
covars <- martini:::arrange_covars(minigwas, covars)

test_that("chi2 are correctly calculated", {
  
  true_vals <- c(98, 98, 97, 99, 97, 96, 0, 1, 1, 1, 1, 0, 
                 0, 0, 99, 99, 97, 98, 98, 99, 0, 1, 0, 1, 0)
  names(true_vals) <- minigwas[['map']][['snp.names']]
  c <- single_snp_association(minigwas, covars, 'chi2')
  
  expect_equal(c, true_vals, tolerance = .1)
  
  # covariates have no effect
  c <- single_snp_association(minigwas, data.frame(), 'chi2')
  expect_equal(c, true_vals, tolerance = .1)
  
})

test_that("glm's chi2 are correctly calculated", {
  
  # without covariates
  true_vals <- c(99, 99, 98, 100, 98, 97, 0, 1, 1, 1, 1, 0, 0, 
                 0, 100, 100, 98, 99, 99, 100, 0.3, 1, 0, 1, 0)
  names(true_vals) <- minigwas[['map']][['snp.names']]
  
  c <- single_snp_association(minigwas, data.frame(), 'glm')
  
  expect_equal(c, true_vals, tolerance = .1)
  
  # with explanatory covariates
  true_vals <- rep(0, nrow(minigwas[['map']]))
  names(true_vals) <- minigwas[['map']][['snp.names']]
  c <- single_snp_association(minigwas, covars, 'glm')
  
  expect_equal(c, true_vals, tolerance = .1)
  
})