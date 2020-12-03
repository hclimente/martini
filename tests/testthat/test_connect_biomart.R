testthat::skip_on_bioc()

test_that("errors if a species is not in ensembl", {
  
  expect_error(martini:::connect_biomart("tsuccinifaciens"), 
               "unable to find an appropriate database for tsuccinifaciens.")
  
})