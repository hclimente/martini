library(martini)
source("big_network.R")

# we make a cluster with A and another with D
cones <- gwas$map
colnames(cones) <- c("chr","snp","cm","pos","allele.1", "allele.2")
cones$c <- 0
cones$selected <- FALSE
cones$selected[grepl("A", cones$snp)] <- TRUE
cones$selected[grepl("D", cones$snp)] <- TRUE

modules <- get_snp_modules(cones, gi)

test_that("output has the right dimensions", {
  expect_equal(dim(modules), dim(cones) + c(0,1))
  expect_true("module" %in% colnames(modules))
})

test_that("snp order is kept", {
  expect_true(all(modules$snp == cones$snp))
})

test_that("we recover the designated modules", {
  expect_equal(unique(modules$module), c(1, NA, 2))
  expect_equal(unique(modules$module[grepl("A", modules$snp)]), 1)
  expect_equal(unique(modules$module[grepl("D", modules$snp)]), 2)
})
