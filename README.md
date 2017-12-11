# martini

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.824643.svg)](https://doi.org/10.5281/zenodo.824643)
[![Build Status](https://travis-ci.org/hclimente/martini.svg?branch=master)](https://travis-ci.org/hclimente/martini)
[![codecov](https://codecov.io/gh/hclimente/martini/branch/master/graph/badge.svg)](https://codecov.io/gh/hclimente/martini)

martini is an R package to perform GWAS experiment that considers prior biological knowledge. This knowledge is modeled as a network of SNPs, were edges represent functional relationships between them (e.g. belonging to the same gene). Then, it looks for regions of the network associated with the phenotype using [SConES](https://academic.oup.com/bioinformatics/article/29/13/i171/198210).

# Installation

To speed up computation, martini uses the [gin](https://github.com/hclimente/gin) C++ library, So, the first step is installing it. Then open an R terminal and install it like any Bioconductor package (requires `devtools`, do `install.packages('devtools')` if you don`t have it):

```
source("https://bioconductor.org/biocLite.R")
biocLite("hclimente/martini")
```

# Usage

After the installation, please refer to the vignettes for usage examples.