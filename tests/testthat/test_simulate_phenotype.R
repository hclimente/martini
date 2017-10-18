library(martini)
data(examplegwas)

set.seed(0)

# get 20 random snps and 20 random effect sizes
causal <- sample(igraph::V(examplegwas$net), 20)
eff <- rnorm(length(causal))
X <- as(examplegwas$gwas$genotypes, "numeric")
X_causal <- X[, examplegwas$gwas$map$snp.names %in% names(causal)]

# qualitative phenotype
Y_ql <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, 
                           effectSize = eff, 
                           qualitative = TRUE, ncases = 1500, ncontrols = 1500)
# quantitative phenotype
Y_qt <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff)
# not a phenotype for everyone
Y_nas <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff,
                            qualitative = T, ncases = 1400, ncontrols = 1350)$fam$affected

test_that("simulate_phenotype runs", {
  expect_type(Y_ql, "list")
  expect_type(Y_qt, "list")
})


test_that("we get the requested number of phenotypes", {
  expect_equal(sum(is.na(Y_ql)), 0)
  expect_equal(sum(is.na(Y_qt)), 0)
  expect_equal(sum(is.na(Y_nas)), 250)
})

qt_p <- apply(X_causal, 2, function(x){
  wilcox.test(Y_qt$fam$affected, x)$p.value
})

ql_p <- apply(X_causal, 2, function(x){
  df <- data.frame(p = Y_ql$fam$affected, g = x)
  chsq <- chisq.test(table(df))
  chsq$p.value
})

test_that("there is an association between phenotype and genotype", {
  expect_gt(sum(ql_p < 0.05), 10)
  expect_gt(sum(qt_p < 0.05), 10)
})