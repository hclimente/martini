skip_if(skip_long)

test_that("plot_ideogram doesn't crash", {
  
  cones <- scones.cv(test_gwas, test_gi)
  plot_ideogram(test_gwas, cones)
  
})
