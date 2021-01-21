test_that("chi2 are correctly calculated", {
  
  true_vals <- c(98, 98, 97, 99, 97, 96, 0, 1, 1, 1, 1, 0, 
                 0, 0, 99, 99, 97, 98, 98, 99, 0, 1, 0, 1, 0)
  names(true_vals) <- test_gwas[['map']][['snp.names']]
  c <- snp_test(test_gwas, covars, 'chi2')
  
  expect_equal(c, true_vals, tolerance = .1)
  
  # covariates have no effect
  c <- snp_test(test_gwas, data.frame(), 'chi2')
  expect_equal(c, true_vals, tolerance = .1)
  
})

test_that("glm's chi2 are correctly calculated", {
  
  # without covariates
  true_vals <- c(99, 99, 98, 100, 98, 97, 0, 1, 1, 1, 1, 0, 0, 
                 0, 100, 100, 98, 99, 99, 100, 0.3, 1, 0, 1, 0)
  names(true_vals) <- test_gwas[['map']][['snp.names']]
  
  c <- snp_test(test_gwas, data.frame(), 'glm', "binomial", "logit")
  
  expect_equal(c, true_vals, tolerance = .1)
  
  # with explanatory covariates
  true_vals <- rep(0, nrow(test_gwas[['map']]))
  names(true_vals) <- test_gwas[['map']][['snp.names']]
  c <- snp_test(test_gwas, covars, 'glm', "binomial", "logit")
  
  expect_equal(c, true_vals, tolerance = .1)
  
})

test_that("linear model is correctly calculated", {
  
  # without covariates
  true_vals <- c(99, 99, 98, 100, 98, 97, 0, 1, 1, 1, 1, 0, 0, 
                 0, 100, 100, 98, 99, 99, 100, 0.3, 1, 0, 1, 0)
  names(true_vals) <- test_gwas[['map']][['snp.names']]
  
  c <- snp_test(test_gwas, data.frame(), 'glm', "Gaussian", "identity")
  
  expect_equal(c, true_vals, tolerance = .1)
  
})
