library(martini)
data(examplegwas)

set.seed(0)

# get 20 random snps and 20 random effect sizes
causal <- sample(igraph::V(examplegwas$net), 50)
eff <- rnorm(length(causal))
X <- as(examplegwas$gwas$genotypes, "numeric")
X_causal <- X[, examplegwas$gwas$map$snp.names %in% names(causal)]

# case-control phenotype
sim_cc <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, 
                             effectSize = eff, 
                             qualitative = TRUE, ncases = 1500, ncontrols = 1500, prevalence = 0.5)
# quantitative phenotype
sim <- simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff)

test_that("output is as expected", {
  
  # right class
  expect_type(sim_cc, "list")
  expect_type(sim, "list")
  
  # only original controls are used
  expect_equal(nrow(sim_cc$fam), 3000)
  expect_equal(dim(sim_cc$genotypes), c(3000, 1800))
  expect_equal(nrow(sim$fam), 3000)
  expect_equal(dim(sim$genotypes), c(3000, 1800))
  
  # genotypes are still a SnpMatrix
  expect_equal(class(sim$genotypes)[1], "SnpMatrix")
  expect_equal(class(sim_cc$genotypes)[1], "SnpMatrix")
  
})

test_that("we get the requested number of phenotypes", {
  
  expect_equal(sum(is.na(sim_cc)), 0)
  expect_equal(sum(is.na(sim)), 0)
  
  # less requested cases and controls than the total number of samples
  expect_equal(sum(is.na(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, 
                                            effectSize = eff, qualitative = T, 
                                            ncases = 1400, ncontrols = 1350, 
                                            prevalence = 0.5)$fam$affected)), 250)
  
})

test_that("there is an association between phenotype and genotype", {
  
  cc_p <- apply(X_causal, 2, function(x){
    df <- data.frame(p = sim_cc$fam$affected, g = x)
    chsq <- chisq.test(table(df))
    chsq$p.value
  })
  expect_gt(sum(cc_p < 0.05), 20)
  
  qt_p <- apply(X_causal, 2, function(x){
    wilcox.test(sim$fam$affected, x)$p.value
  })
  expect_gt(sum(qt_p < 0.05), 20)
  
})

test_that("errors when it should", {
  
  source("big_network.R")
  badCausal <- sample(igraph::V(gi), 10)
  
  # general errors
  expect_error(simulate_phenotype(examplegwas$gwas, badCausal, h2 = 1, effectSize = eff),
               paste("The following causal SNPs are not in the SNP list:", 
                     paste(names(badCausal), collapse = ",")), fixed=TRUE)
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = 1.5, effectSize = eff),
               "h2 must be between 0 and 1. Current value is 1.5.", fixed=TRUE)
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = -1.5, effectSize = eff),
               "h2 must be between 0 and 1. Current value is -1.5.", fixed=TRUE)
  
  # errors related to qualitative phenotypes
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff, qualitative = T), 
               'argument "ncases" is missing, with no default', fixed=TRUE)
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff, qualitative = T, ncases = 1e4), 
               'argument "ncontrols" is missing, with no default', fixed=TRUE)
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff, qualitative = T, ncases = 1500, ncontrols = 1500), 
               'argument "prevalence" is missing, with no default', fixed=TRUE)
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff, qualitative = T, ncases = 1e4, ncontrols = 1e4, prevalence = 0.1), 
               "Requested number of cases and controls too high (> # samples).", fixed=TRUE)
  expect_error(simulate_phenotype(examplegwas$gwas, causal, h2 = 1, effectSize = eff, qualitative = T, ncases = 1500, ncontrols = 1500, prevalence = 0.1), 
               "Requested number of cases too high (> # samples * prevalence).", fixed=TRUE)
  
})