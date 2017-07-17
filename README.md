# martini

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.824643.svg)](https://doi.org/10.5281/zenodo.824643)
[![Build Status](https://travis-ci.org/hclimente/martini.svg?branch=master)](https://travis-ci.org/hclimente/martini)

martini is an R interface of [gin](https://github.com/hclimente/gin). gin performs GWAS incorporating prior knowledge, namely biological networks. martini provides an R interface for gin's signature `shake` function, and extends it useful functions for data preprocessing, and plotting and analyzing the results.

# Installation

If you haven't, first, install [`gin`](https://github.com/hclimente/gin). Then open an R terminal and install it like any Bioconductor package (requires `devtools`, do `install.packages('devtools')` if you don`t have it):

```
source("https://bioconductor.org/biocLite.R")
biocLite("hclimente/martini")
```

# Usage

## GWAS Incorpotating Networks

```{r}
library(martini)
data(examplegwas)
```
The example data contains two variables:

- `gwas`: a `snpMatrix` structure, containing the genotype information (e.g. created by `read.pedfile()`, from the `snpStats` package).
- `net`: an `igraph` network, containing the SNP-network (e.g. created by `get_GS_network()`).

martini uses igraph networks of SNPs. The user can connect them according to gene membership, sequence contiguity, protein-protein interactions, etc.

```{r}
g <- shake(gwas, net)
```

`shake` is the main function in martini. Additional arguments can be passed. `shake` returns a copy of the `map` from the `gwas` object, which contain information of the SNPs (name, chromosome, genomic position...). `shake` adds three columns.

- `C` is a numeric vector with the association scone for each SNP.
- `selected` is a boolean vector informing about if the SNP was selected or not by `shake`.
- `cluster` is a integer vector with information about which SNPs is adjacent to which SNP in the network. The integer is the identifier of that cluster (NA if the SNP is disconnected from any other selected SNP).

```{r}
head(g)
#   V1 snp.names V3 V4 allele.1 allele.2         C  selected  cluster
# 1  1       1_1  0  1        A        T 361.13735     FALSE       NA
# 2  1       1_2  0  2        T        A 344.29586     FALSE       NA
# 3  1       1_3  0  3        T        A 894.68186     FALSE        1
# 4  1       1_4  0  4        T        A 180.98245     FALSE        1
# 5  1       1_5  0  5        A        T 402.55416     FALSE        1
# 6  1       1_6  0  6        A        T  21.54136     FALSE        1

```

## Simulate quantitative phenotype

```{r}
library(martini)
data(examplegwas)

simulation <- data.frame(causal = logical(1800) )

# simulate 20 causal SNPs, interconnected in the PPI network
simulation$causal <- simulate_causal_snps(gwas, net, 20)

# get their effect sizes from a normal distribution and simulate the phenotype
simulation$effectSize[simulation$causal] <- rnorm(sum(simulation$causal))
Y.simu <- simulate_phenotype(gwas, simulation$causal, 
                            h2 = 1, 
                            effectSize = simulation$effectSize[simulation$causal], 
                            qualitative = TRUE, ncases = 3000, ncontrols = 3000)

# study the association between the SNPs and the phenotype
simulation$pval <- apply(as(gwas$genotypes, "numeric"), 2, function(x){
    df <- data.frame(p = Y.simu, g = x)
    chsq <- chisq.test(table(df))
    chsq$p.value
  })

ggplot(simulation, aes(pval, fill = causal)) +
  geom_histogram()

ggplot(subset(simulation, causal), aes(x = abs(effectSize), y = -log10(pval))) +
  geom_point() +
  geom_smooth(method='lm')

```
