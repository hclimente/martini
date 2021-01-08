bed <- gwas2bed(test_gwas)

test_that("output is as expected", {
  
  groups <- group_snps(bed, 'chr', 'start', 5)
  
  expect_type(groups, 'list')
  expect_type(groups[,1], 'character')
  expect_type(groups[,2], 'integer')
  expect_type(groups[,3], 'integer')
  expect_equal(ncol(groups), 3)
  
})

test_that("we do group snps", {
  
  groups <- group_snps(bed, 'chr', 'start', 5)
  expect_true(any(groups[,1] == 'chr1' & groups[,2] == 60 & groups[,3] == 74))
  expect_true(any(groups[,1] == 'chr1' & groups[,2] == 110 & groups[,3] == 124))
  expect_true(any(groups[,1] == 'chr2' & groups[,2] == 35 & groups[,3] == 39))
  expect_true(any(groups[,1] == 'chr2' & groups[,2] == 75 & groups[,3] == 89))
  
})

test_that("threshold works properly", {
  
  # all SNPs are at distance <= 10 of another one
  groups <- group_snps(bed, 'chr', 'start', 11)
  expect_equal(nrow(unique(groups)), 2)
  
  # some SNPs are at distance 5 of another one
  groups <- group_snps(bed, 'chr', 'start', 5)
  expect_equal(nrow(unique(groups)), nrow(bed) - 6)
  
  # differences dissapear when going below 5
  groups <- group_snps(bed, 'chr', 'start', 4)
  expect_equal(nrow(unique(groups)), nrow(bed))
  
  # all SNPs are separated
  groups <- group_snps(bed, 'chr', 'start', 1)
  expect_equal(nrow(unique(groups)), nrow(bed))
  
})
