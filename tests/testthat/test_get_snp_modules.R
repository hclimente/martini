cones <- scones.cv(test_gwas, test_gm)
modules <- get_snp_modules(test_gwas, cones)

test_that("output has the right dimensions", {
  expect_equal(dim(modules), dim(test_gwas[['map']]) + c(0,1))
  expect_true("module" %in% colnames(modules))
})

test_that("snp order is kept", {
  expect_true(all(modules$snp == cones$snp))
})

test_that("we recover the designated modules", {
  expect_equal(unique(modules$module), c(1, NA, 2))
  expect_equal(unique(modules$module[grepl("A", modules$snp)]), 1)
  expect_equal(unique(modules$module[grepl("C", modules$snp)]), 2)
})
