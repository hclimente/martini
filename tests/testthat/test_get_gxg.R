testthat::skip_on_bioc()
library(martini)

good_dog <- martini:::get_gxg_biogrid(9615)
good_frog <- martini:::get_gxg_string(8364)

start <- proc.time()
suppressWarnings(dog <- martini:::get_gxg('biogrid', 9615, flush = TRUE))
t_dog <- proc.time() - start

start <- proc.time()
suppressWarnings(frog <- martini:::get_gxg('string', 8364, flush = TRUE))
t_frog <- proc.time() - start

test_that("output is as expected", {
  
  expect_equal(good_dog, dog)
  expect_equal(good_frog, frog)
  
})

test_that("cache works", {

  start <- proc.time()
  suppressWarnings(dog_cached <- martini:::get_gxg('biogrid', 9615, flush = FALSE))
  t_dog_cached <- proc.time() - start
  
  start <- proc.time()
  suppressWarnings(frog_cached <- martini:::get_gxg('string', 8364, flush = FALSE))
  t_frog_cached <- proc.time() - start
  
  # results are the same
  expect_equal(dog, dog_cached)
  expect_equal(frog, frog_cached)
  
  # elapsed time is reduced when we cache results
  expect_gt(t_dog['elapsed'], t_dog_cached['elapsed'])
  expect_gt(t_frog['elapsed'], t_frog_cached['elapsed'])
  
  # appropriate warnings
  expect_warning(martini:::get_gxg('biogrid', 9615, flush = FALSE),
                 "using cache. Use flush = TRUE to get new gene interactions.")
  expect_message(martini:::get_gxg('biogrid', 9615, flush = TRUE),
                 "cache flushed!")

})