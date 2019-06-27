# martini

[![Build Status](https://travis-ci.org/hclimente/martini.svg?branch=master)](https://travis-ci.org/hclimente/martini)
[![Build status](https://ci.appveyor.com/api/projects/status/ehnhhg2n5cs7pkk6?svg=true)](https://ci.appveyor.com/project/hclimente/martini)
[![codecov](https://codecov.io/gh/hclimente/martini/branch/master/graph/badge.svg)](https://codecov.io/gh/hclimente/martini)
[![BioC](https://bioconductor.org/shields/years-in-bioc/martini.svg)](https://bioconductor.org/packages/devel/bioc/html/martini.html)

`martini` is an R package to perform GWAS experiment that considers prior biological knowledge. This knowledge is modeled as a network of SNPs, were edges represent functional relationships between them (e.g. belonging to the same gene). Then, it looks for regions of the network associated with the phenotype using [SConES](https://academic.oup.com/bioinformatics/article/29/13/i171/198210).

# Installation

Install `martini` like any Bioconductor package:

``` r
install.packages("BiocManager")
BiocManager::install("martini")
```

# Usage

Running `martini` is a three step process:

``` r
library(martini)

# 1. Read GWAS data with read.pedfile (or load the example :) )
data(minigwas)

# 2. Create the SNP network: GS (structural information), GM (GS + gene 
# annotation information) or GI (GM + protein-protein interaction information)
gs <- get_GS_network(minigwas)

# 3. Find connected, explanatory SNPs (cones)
cones <- search_cones(minigwas, gs)

# cones$selected informs about whether the SNP is selected as cones or not
head(cones)
#   snp chr cm pos allele.1 allele.2        c selected module
# 3 1A1   1  0  10        A        G 96.15385     TRUE      1
# 4 1A2   1  0  20        A        G 96.15385     TRUE      1
# 5 1A3   1  0  30        A        G 96.15385     TRUE      1
# 6 1A4   1  0  40        A        G 96.15385     TRUE      1
# 7 1A5   1  0  50        A        G 96.15385     TRUE      1
# 8 1A6   1  0  60        A        G 96.15385     TRUE      1
```

Please, refer to the vignettes for more detailed usage examples. `martini` results can be further examined using the [blur](https://github.com/hclimente/blur) package.
