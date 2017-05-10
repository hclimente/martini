# rscones

R interface for SConES.

## Installation

Currently only local installation is supported. First, clone the repository with all its submodules.

```{bash}
git clone --recursive git@github.com:hclimente/rscones.git
cd rscones/src/scones2
git submodule update --init --recursive
cd ../..
```

Then, open a R terminal and install it from source:

```{r}
install.packages(".", repos = NULL, type="source")
```
