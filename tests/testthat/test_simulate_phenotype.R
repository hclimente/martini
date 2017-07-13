library(martini)
data(examplegwas)

set.seed(0)

test_that("simulate_phenotype runs", {
  expect_type(simulate_phenotype(gwas, c("rs1", "rs2"), h2 = 1), "double")
  expect_type(simulate_phenotype(gwas, c("rs1", "rs2"), h2 = 1, qualitative = T, ncases = 3000, ncontrols = 3000), "double")
})