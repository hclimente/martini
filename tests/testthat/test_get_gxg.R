library(martini)

# retrieve PPI from the dog, which in Biogrid 3.4.151 had 25 described interactions
good_dog <- martini:::get_gxg_biogrid(9615)

start <- proc.time()
dog <- martini:::get_gxg('biogrid', 9615, flush = TRUE)
elapsed <- proc.time() - start

test_that("output is as expected", {
  
  # output
  expect_equal(good_dog, dog)
  
})

test_that("cache works", {

  start <- proc.time()
  dog_cached <- martini:::get_gxg('biogrid', 9615, flush = FALSE)
  elapsed_cached <- proc.time() - start
  
  # results are the same
  expect_equal(dog, dog_cached)
  
  # elapsed time is reduced when we cache results
  expect_gt(elapsed['elapsed'], elapsed_cached['elapsed'])
  
  # appropriate warnings
  expect_warning(martini:::get_gxg('biogrid', 9615, flush = FALSE),
                 "using cache. Use flush = TRUE to get new gene interactions.")
  expect_message(martini:::get_gxg('biogrid', 9615, flush = TRUE),
                 "cache flushed!")

})