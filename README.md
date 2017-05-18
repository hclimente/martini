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
