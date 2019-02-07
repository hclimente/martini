library(martini)

set.seed(0)

test_that("default args are being set", {
  expect_equal(parse_scones_settings(c = 1)$criterion, 'consistency')
  expect_equal(parse_scones_settings(c = 1)$score, 'chi2')
  expect_equal(parse_scones_settings(c = 1)$etas, rep(1, 10))
  expect_equal(parse_scones_settings(c = 1)$lambdas, rep(1, 10))
})

test_that("default we can change values", {
  expect_equal(parse_scones_settings(c = 1, criterion = "consistency")$criterion, "consistency")
  expect_equal(parse_scones_settings(c = 1, criterion = "bic")$criterion, 'bic')
  expect_equal(parse_scones_settings(c = 1, criterion = "aic")$criterion, 'aic')
  expect_equal(parse_scones_settings(c = 1, criterion = "aicc")$criterion, 'aicc')
  expect_error(parse_scones_settings(c = 1, criterion = "kk"))
  expect_equal(parse_scones_settings(c = 1, score = "glm")$score, 'glm')
  expect_equal(parse_scones_settings(c = 1, score = "chi2")$score, 'chi2')
  expect_error(parse_scones_settings(c = 1, score = "kk"))
  expect_equal(parse_scones_settings(c = 1, etas = c(3,4,5))$etas, c(3,4,5))
  expect_error(parse_scones_settings(etas = c("a","b"))$etas)
  expect_equal(parse_scones_settings(c = 1, lambdas = c(3,4,5))$lambdas, c(3,4,5))
  expect_error(parse_scones_settings(lambdas = c("a","b"))$lambdas)
  expect_error(parse_scones_settings(c = 1, debug = 3))
})