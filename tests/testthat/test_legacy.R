test_that("search_cones calls the right function", {
 
  expect_warning(cones <- search_cones(test_gwas, test_gi),
                 "'search_cones' is deprecated.")
  expect_equal(cones, scones.cv(test_gwas, test_gi))
  expect_equal(suppressWarnings(search_cones(test_gwas, test_gi, sigmod = TRUE)), 
               sigmod.cv(test_gwas, test_gi))
  expect_warning(search_cones(test_gwas, test_gi, encoding = 'foo'),
                 "Argument encoding ignored, will be set to additive.")
  
})
