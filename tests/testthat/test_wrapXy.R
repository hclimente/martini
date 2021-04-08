test_that("output is as expected", {
  
  i <- martini:::wrap_Xy(X, Y, test_gwas$map$snp.name, test_gi)
  
  expect_equal(i[['gwas']][['genotypes']], X)
  expect_equal(i[['gwas']][['fam']]$affected, Y)
  expect_equal(i[['gwas']][['map']]$snp, test_gwas$map$snp.name)
  
  expect_true(all(V(test_gi) == V(i[['net']])))
  expect_true(all(E(test_gi) == E(i[['net']])))
  
})
