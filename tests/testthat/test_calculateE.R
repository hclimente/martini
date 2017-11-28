library(martini)

test_that("E has the expected output", {
  
  expect_equal(length(calculateE(rnorm(322), 0.25)), 322)
  expect_equal(length(calculateE(rnorm(7), 0.25)), 7)
  expect_equal(class(calculateE(rnorm(7), 0.25)), "numeric")
  
})

test_that("E depends on variance in genotypes", {
  
  expect_equal(calculateE(rep(1,10), 1), rep(0,10))
  expect_equal(calculateE(rep(5.5,10), 1), rep(0,10))
  
})

test_that("E depends on heritability", {
  
  G <- rnorm(100)
  expect_true(var(calculateE(G, 0.25)) > var(calculateE(G, 0.75)))
  
})