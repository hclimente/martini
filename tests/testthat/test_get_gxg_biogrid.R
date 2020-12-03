library(martini)

# retrieve PPI from the dog, which in Biogrid 3.4.151 had 25 described interactions
dog <- martini:::get_gxg_biogrid(9615)

test_that("output is as expected", {
  
  # output
  expect_true(is.data.frame(dog))
  expect_equal(ncol(dog), 2)
  expect_equal(colnames(dog), c("gene1","gene2"))
  
  # we get the symbols
  expect_true("NEDD4" %in% c(dog$gene1, dog$gene2))
  expect_true("CREBBP" %in% c(dog$gene1, dog$gene2))
  expect_true("SP1" %in% c(dog$gene1, dog$gene2))
  
})

test_that("we retrieve known interactions", {
  
  # known interactions
  expect_equal(sum(dog$gene1 == "CREBBP" & dog$gene2 == "SP1"), 1)
  expect_equal(sum(dog$gene1 == "DTNA" & dog$gene2 == "ACTB"), 1)
  
  # expected order
  expect_equal(sum(dog$gene1 == "SP1" & dog$gene2 == "CREBBP"), 0)
  expect_equal(sum(dog$gene1 == "ACTB" & dog$gene2 == "DTNA"), 0)
  
})