test_that("subsets the right columns", {
 
  test_snpMapping <- data.frame(a = c('geneA'), b = c('trash'), c = c('snp1'))
  sanitized <- sanitize_snpMapping(test_snpMapping, c('c','a'))
  
  expect_equal(c(1,2), dim(sanitized))
  expect_equal(c('snp','gene'), colnames(sanitized))
  expect_equal('geneA', sanitized[1,'gene'])
  expect_equal('snp1', sanitized[1,'snp'])
  
})

test_that("warnings are correct", {
  
  expect_warning(sanitize_snpMapping(NULL, NULL),
                 "no mappings between SNPs and genes were provided.")
  expect_equal(c(0,2), dim(suppressWarnings(sanitize_snpMapping(NULL, NULL))))
  
})