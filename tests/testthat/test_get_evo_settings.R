library(martini)
data(examplegwas)

set.seed(0)

test_that("default args are being set", {
  expect_equal(get_evo_settings()$encoding, 0)
  expect_equal(get_evo_settings()$modelScore, 1)
  expect_equal(get_evo_settings()$associationScore, 0)
})

test_that("default we can change values", {
  expect_equal(get_evo_settings(encoding = "additive")$encoding, 0)
  expect_equal(get_evo_settings(encoding = "recessive")$encoding, 1)
  expect_equal(get_evo_settings(encoding = "dominant")$encoding, 2)
  expect_equal(get_evo_settings(encoding = "codominant")$encoding, 3)
  expect_equal(get_evo_settings(modelScore = "consistency")$modelScore, 0)
  expect_equal(get_evo_settings(modelScore = "bic")$modelScore, 1)
  expect_equal(get_evo_settings(modelScore = "aic")$modelScore, 2)
  expect_equal(get_evo_settings(modelScore = "aicc")$modelScore, 3)
  expect_equal(get_evo_settings(modelScore = "mbic")$modelScore, 4)
  expect_error(get_evo_settings(modelScore = "kk"))
  expect_equal(get_evo_settings(associationScore = "skat")$associationScore, 0)
  expect_equal(get_evo_settings(associationScore = "chi2")$associationScore, 1)
  expect_equal(get_evo_settings(associationScore = "trend")$associationScore, 2)
  expect_error(get_evo_settings(associationScore = "kk"))
})