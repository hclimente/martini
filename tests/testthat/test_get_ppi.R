library(martini)

# retrieve PPI from the dog, which in Biogrid 3.4.151 had 25 described interactions
dog <- get_ppi(9615)

test_that("output is as expected", {
  
  # output
  expect_true(is.data.frame(dog))
  expect_equal(ncol(dog), 2)
  expect_equal(colnames(dog), c("geneA","geneB"))
  
  # we get the symbols
  expect_true("NEDD4" %in% c(dog$geneA, dog$geneB))
  expect_true("CREBBP" %in% c(dog$geneA, dog$geneB))
  expect_true("SP1" %in% c(dog$geneA, dog$geneB))

})

test_that("we retrieve known interactions", {

  # known interactions
  expect_equal(sum(dog$geneA == "CREBBP"	& dog$geneB == "SP1"), 1)
  expect_equal(sum(dog$geneA == "DTNA"	& dog$geneB == "ACTB"), 1)
  
  # expected order
  expect_equal(sum(dog$geneA == "SP1"	& dog$geneB == "CREBBP"), 0)
  expect_equal(sum(dog$geneA == "ACTB"	& dog$geneB == "DTNA"), 0)

})