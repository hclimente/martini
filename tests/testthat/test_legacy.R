test_that("output is as expected", {
  
  suppressWarnings(cones <- search_cones(test_gwas, test_gi))
  
  expect_equal(dim(cones), dim(test_gwas$map) + c(0,3))
  expect_equal(class(cones), "data.frame")
  expect_equal(class(cones$selected), "logical")
  expect_equal(class(cones$c), "numeric")
  expect_equal(class(cones$module), "numeric")
  
})

test_that("search_cones calls the right function", {
 
  expect_warning(cones <- search_cones(test_gwas, test_gi),
                 "'search_cones' is deprecated.")
  expect_equal(sort(cones$snp[cones$selected]),
               sort(names(V(scones.cv(test_gwas, test_gi)))))
  
  suppressWarnings(cones <- search_cones(test_gwas, test_gi, sigmod = TRUE))
  expect_equal(sort(cones$snp[cones$selected]),
               sort(names(V(sigmod.cv(test_gwas, test_gi)))))
  
  expect_warning(search_cones(test_gwas, test_gi, encoding = 'foo'),
                 "Argument encoding ignored, will be set to additive.")
  
})
