test_that("output is as expected", {
  
  cones <- sigmod.cv(test_gwas, test_gi)
  
  expect_equal(class(cones), "igraph")
  expect_true(all(names(V(cones)) %in% test_map$snp))
  
})

test_that("we recover causal SNPs", {
  
  cones <- sigmod.cv(test_gwas, test_gi, etas = 0, lambdas = 0)
  
  # wrong eta and lambda return the trivial solution
  expect_equal(gorder(cones), gorder(test_gi))
  
  set.seed(42)
  cones <- sigmod.cv(test_gwas, test_gi, 
                     etas = seq(2, 0, length=10), 
                     lambdas = seq(2, 0, length=10))
  
  skip_on_os("windows")
  expect_equal(gorder(cones), sum(grepl("[AC]", test_map$snp)))
  
  set.seed(42)
  cones <- sigmod.cv(test_gwas, test_gi, 
                     etas = seq(2, 0, length=10), 
                     lambdas = seq(2, 0, length=10),
                     criterion = 'bic')
  
  skip_on_os("windows")
  expect_equal(gorder(cones), sum(grepl("[AC]", test_map$snp)))
  expect_equal(sort(names(V(cones))), 
               sort(grep("[AC]", test_map$snp, value = TRUE)))
  
})
