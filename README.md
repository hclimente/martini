# martini

R implementation of [gin](https://github.com/hclimente/gin).

# Installation

Currently only local installation is supported. First, clone the repository with all its submodules.

```{bash}
git clone --recursive git@github.com:hclimente/martini.git
cd martini/src/gin
git submodule update --init --recursive
cd ../..
```

Then, open a R terminal and install it from source:

```{r}
install.packages(".", repos = NULL, type="source")
```

# Usage

## GWAS Incorpotating Networks

```{r}
library(martini)

data(simplegwas)
# contains two variables:
#   - geno: genotype information, created with read.pedfile (snpStats package)
#   - net: dataframe with edge information

graph <- graph_from_data_frame(net, directed = F )

g <- shake(geno, graph)

# shake creates a copy from the geno object with two modifications
# g$map is modified to include two extra columns: one with the statistic and other informing if the feature got selected

head(g$map)
#   V1 snp.names V3 V4 allele.1 allele.2  ginscore ginpicked
# 1  1       1_1  0  1        A        T 361.13735     FALSE
# 2  1       1_2  0  2        T        A 344.29586     FALSE
# 3  1       1_3  0  3        T        A 894.68186     FALSE
# 4  1       1_4  0  4        T        A 180.98245     FALSE
# 5  1       1_5  0  5        A        T 402.55416     FALSE
# 6  1       1_6  0  6        A        T  21.54136     FALSE

# g$gin contains the best parameters found in the gridsearch

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
