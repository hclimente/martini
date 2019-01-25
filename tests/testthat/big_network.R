library(igraph)

# MAKE GWAS OBJECT
## map
gwas <- list()
gwas$map <- read.table(text = "
                       chr snp.names cm gpos allele.1 allele.2
                       1 1A1 0 10 A G
                       1 1A2 0 20 A G
                       1 1A3 0 30 A G
                       1 1A4 0 40 A G
                       1 1A5 0 50 A G
                       1 1A6 0 60 A G
                       1 1-1 0 65 A G
                       1 1B1 0 70 A G
                       1 1B2 0 80 A G
                       1 1B3 0 90 A G
                       1 1B4 0 100 A G
                       1 1B5 0 110 A G
                       1 1B6 0 115 A G
                       1 1-2 0 120 A G
                       2 2C1 0 35 A G
                       2 2C2 0 45 A G
                       2 2C3 0 55 A G
                       2 2C4 0 65 A G
                       2 2C5 0 75 A G
                       2 2C6 0 80 A G
                       2 2-1 0 85 A G
                       2 2-2 0 95 A G
                       2 2D1 0 105 A G
                       2 2D2 0 115 A G
                       2 2D3 0 125 A G
                       ", header = TRUE, stringsAsFactors = FALSE)

## genotypes
N <- 100
sol <- grepl("[AC]", gwas$map$snp.names)
pCausal <- sum(sol)
pNonCausal <- nrow(gwas$map) - pCausal
causal <- c(rep(2, N/2), rep(0, N/2))
rest <- rep(0,N)

X <- do.call(cbind, lapply(sol, function(x) if(x) causal else rest))
colnames(X) <- gwas$map$snp.names
rownames(X) <- 1:nrow(X)
Xp <- X + 1
mode(Xp) <- "raw"
gwas$genotypes <- new("SnpMatrix", Xp)

## phenotypes
Y <- c(rep(2, N/2), rep(1, N/2))
gwas$fam <- data.frame(pedigree = 1:N,
                       member = 1:N,
                       father = NA,
                       mother = NA,
                       sex = sample(c(1,2), N, replace = TRUE),
                       affected = Y)

## MAKE NETWORKS
snpMapping <- data.frame(snp = gwas$map$snp.names,
                         gene = substr(gwas$map$snp.names, 2, 2))
snpMapping <- subset(snpMapping, gene != "-")
ppi <- read.table(text = "
                  gene1 gene2
                  A B
                  A C
                  B D
                  ", header = TRUE, stringsAsFactors = FALSE)

gs <- get_GS_network(gwas)
gm <- get_GM_network(gwas, snpMapping = snpMapping)
gi <- get_GI_network(gwas, snpMapping = snpMapping, ppi = ppi)
W <- as_adj(gi)

solution <- search_cones(gwas, gi)

# CREATE MINI
# minigwas = gwas
# minigs = gs
# minigm = gm
# minigi = gi
# miniX = X
# miniY = Y
# miniW = W
# minisnpMapping = snpMapping
# minippi = ppi
# miniSolution = solution
# 
# save(minigwas, file = "data/minigwas.rda")
# save(minigs, file = "data/minigs.rda")
# save(minigm, file = "data/minigm.rda")
# save(minigi, file = "data/minigi.rda")
# save(miniX, file = "data/miniX.rda")
# save(miniY, file = "data/miniY.rda")
# save(miniW, file = "data/miniW.rda")
# save(minisnpMapping, file = "data/minisnpMapping.rda")
# save(minippi, file = "data/minippi.rda")
# save(miniSolution, file = "data/miniSolution.rda")