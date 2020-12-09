library(martini)

test_that("detects when a package is not installed", {
  
  expect_error(check_installed("fake package 123"), "fake package 123 needed for this function to work. Please install it.")
  expect_error(check_installed("fake package 123", "fake fun"), "fake package 123 needed for fake fun to work. Please install it.")
  expect_error(check_installed("martini"), NA)
  
})
