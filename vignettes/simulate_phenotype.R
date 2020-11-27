## ----load_martini-------------------------------------------------------------
library(martini)

## ----create_gi----------------------------------------------------------------
gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)

## ----simulate_snps------------------------------------------------------------
causal <- simulate_causal_snps(gi, ngenes = 2, pcausal = 0.3)

par(mar=c(0,0,0,0)+.1)
plot(gi, mark.groups = names(causal))

## ----simulate_quantitative_phenotype------------------------------------------
# create a random 100x25 matrix of genotypes
X <- lapply(seq_len(100), function(i) { sample(c(0,1,2), 25, replace = TRUE) })
X <- do.call(rbind, X)

colnames(X) <- minigwas$map$snp.names
rownames(X) <- 1:nrow(X)

X <- X + 1
mode(X) <- "raw"
minigwas$genotypes <- new("SnpMatrix", X)

simulated <- simulate_phenotype(minigwas, snps = causal, h2 = 0.9)

## ----simulate_qualitative_phenotype-------------------------------------------
simulated <- simulate_phenotype(minigwas, snps = causal, h2 = 0.9, qualitative = TRUE,
                                ncases = 10, ncontrols = 40, prevalence = 0.2)

