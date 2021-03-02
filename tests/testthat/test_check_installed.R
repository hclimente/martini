library(martini)

test_that("detects when a package is not installed", {
  
  expect_error(check_installed(c("foo")), "This function requires the following packages to be installed:\nfoo")
  expect_error(check_installed(c("foo","bar"), "baz"), "baz requires the following packages to be installed:\nfoo\nbar")
  expect_error(check_installed("martini"), NA)
  
})
