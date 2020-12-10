library(martini)

gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)

test_that("output is as expected", {
  
  expect_equal(sigmod(minigwas, gi, 2.25, 0.034), sigmod.cv(minigwas, gi))
  
})
