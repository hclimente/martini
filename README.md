# martini

martini is an R interface of [gin](https://github.com/hclimente/gin). gin performs GWAS incorporating prior knowledge, namely biological networks. martini provides an R interface for gin's signature `shake` function, and extends it useful functions for data preprocessing, and plotting and analyzing the results.

# Installation

First, install [gin](https://github.com/hclimente/gin). Then,

```
devtools::install_github("hclimente/martini")
```

# Usage

## GWAS Incorpotating Networks

```{r}
library(martini)
data(simplegwas)
```
The example data contains two variables:

- `geno`: genotype information, created with read.pedfile (`snpStats` package)
- `net`: dataframe with edge information.

```{r}
graph <- graph_from_data_frame(net, directed = F )
```

martini uses igraph networks of SNPs. The user can connect them according to gene membership, sequence contiguity, protein-protein interactions, etc.

```{r}
g <- shake(geno, graph)
```

`shake` is the main function in martini. Additional arguments can be passed. `shake` creates a copy from the `geno` object with two modifications:

- `g$map` is modified to include two extra columns: one with the statistic and other informing if the feature got selected.
- `g$gin` contains the best gin parameters found in the gridsearch.

```{r}
head(g$map)
#   V1 snp.names V3 V4 allele.1 allele.2  ginscore ginpicked
# 1  1       1_1  0  1        A        T 361.13735     FALSE
# 2  1       1_2  0  2        T        A 344.29586     FALSE
# 3  1       1_3  0  3        T        A 894.68186     FALSE
# 4  1       1_4  0  4        T        A 180.98245     FALSE
# 5  1       1_5  0  5        A        T 402.55416     FALSE
# 6  1       1_6  0  6        A        T  21.54136     FALSE

g$gin
# $lambda
# [1] 278.2559
# 
# $eta
# [1] 16681.01

```

## Simulate quantitative phenotype

```{r}
library(martini)
data(toyGWAS)

simulation <- data.frame(causal = logical(1000) )

# simulate 20 causal SNPs, interconnected in the PPI network
simulation$causal <- simulateCausalSNPs(toyGWAS$net, 20)

# get their effect sizes from a normal distribution and simulate the phenotype
simulation$effectSize[simulation$causal] <- rnorm(sum(simulation$causal))
Y.simu <- simulatePhenotype(toyGWAS$X, simulation$causal, 
                            h2 = 1, 
                            effectSize = simulation$effectSize[simulation$causal], 
                            qualitative = TRUE, ncases = 250, ncontrols = 250)

# study the association between the SNPs and the phenotype
simulation$pval <- apply(toyGWAS$X, 2, function(x){
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
