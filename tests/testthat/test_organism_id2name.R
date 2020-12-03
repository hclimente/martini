test_that("we get the right species", {
  expect_equal(martini:::organism_id2name(9606), "hsapiens")
  expect_equal(martini:::organism_id2name(3702), "athaliana")
  expect_equal(martini:::organism_id2name(167), "tsuccinifaciens")
})