library(martini)

test_that("conversions are right", {
  
  additive <- matrix(c(0,1,2,1,2,0,2,1,0), 3, 3)
  recessive <- matrix(c(0,0,1,0,1,0,1,0,0), 3, 3)
  dominant <- matrix(c(0,1,1,1,1,0,1,1,0), 3, 3)
  codominant <- matrix(c(0,1,0,1,0,0,0,1,0), 3, 3)
  
  expect_equal(encode_gwas(additive, "additive"), additive)
  expect_equal(encode_gwas(additive, "recessive"), recessive)
  expect_equal(encode_gwas(additive, "dominant"), dominant)
  expect_equal(encode_gwas(additive, "codominant"), codominant)

})

test_that("output is as expected", {
  
  data(examplegwas)
  Xa <- as(examplegwas$gwas$genotypes, "numeric")
  Xr <- encode_gwas(Xa, "recessive")
  
  expect_equal(dim(Xa), dim(Xr))
  expect_equal(Xa, encode_gwas(Xa, "additive"))
  
})

test_that("errors when it should", {
  
  expect_error(encode_gwas(matrix(), "xyz"), "Invalid encoding.")
  
})