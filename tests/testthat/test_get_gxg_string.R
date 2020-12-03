skip_if(skip_long)

ppi <- martini:::get_gxg_string(9606)

test_that("output is as expected", {
  
  # output
  expect_true(is.data.frame(ppi))
  expect_equal(ncol(ppi), 2)
  expect_equal(colnames(ppi), c("gene1","gene2"))
  
  # we get some known genes
  expect_true("TP53" %in% c(ppi$gene1, ppi$gene2))
  expect_true("AKT1" %in% c(ppi$gene1, ppi$gene2))
  expect_true("GAPDH" %in% c(ppi$gene1, ppi$gene2))
  
})

test_that("we retrieve known interactions", {
  
  # known interactions
  expect_equal(sum(ppi$gene1 == "KLF2" & ppi$gene2 == "TP53"), 1)
  expect_equal(sum(ppi$gene1 == "DUSP6" & ppi$gene2 == "AKT1"), 1)
  
  # expected order
  expect_equal(sum(ppi$gene1 == "TP53" & ppi$gene2 == "KLF2"), 0)
  expect_equal(sum(ppi$gene1 == "AKT1" & ppi$gene2 == "DUSP6"), 0)
  
})