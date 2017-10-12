library(martini)
data(examplegwas)

set.seed(0)
causal <- sample(igraph::V(examplegwas$net), 5)
eff <- rnorm(5)

test_that("simulate_phenotype runs", {
  expect_type(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff), "list")
  expect_type(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff,
                                 qualitative = T, ncases = 1500, ncontrols = 1500), "list")
})

# get 20 random snps and 20 random effect sizes
causal <- sample(igraph::V(examplegwas$net), 20)
eff <- rnorm(length(causal))
X <- as(examplegwas$gwas$genotypes, "numeric")
X_causal <- X[, examplegwas$gwas$map$snp.names %in% names(causal)]

# association between causal SNPs and the qualitative phenotype
Y_ql <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, 
                           effectSize = eff, 
                           qualitative = TRUE, ncases = 1500, ncontrols = 1500)
ql_p <- apply(X_causal, 2, function(x){
  df <- data.frame(p = Y_ql$fam$affected, g = x)
  chsq <- chisq.test(table(df))
  chsq$p.value
})

# association between causal SNPs and the quantitative phenotype
Y_qt <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff)
qt_p <- apply(X_causal, 2, function(x){
  wilcox.test(Y_qt$fam$affected, x)$p.value
})

test_that("there is an association between phenotype and genotype", {
  expect_gt(sum(ql_p < 0.05), 10)
  expect_gt(sum(qt_p < 0.05), 10)
})