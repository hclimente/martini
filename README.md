# martini

![R-CMD-check-bioc](https://github.com/hclimente/martini/workflows/R-CMD-check-bioc/badge.svg)
[![codecov](https://codecov.io/gh/hclimente/martini/branch/master/graph/badge.svg)](https://codecov.io/gh/hclimente/martini)
[![BioC](https://bioconductor.org/shields/years-in-bioc/martini.svg)](https://bioconductor.org/packages/devel/bioc/html/martini.html)

`martini` is an R package to perform GWAS experiment that considers prior biological knowledge. This knowledge is modeled as a network of SNPs, were edges represent functional relationships between them (e.g. belonging to the same gene). Then, it looks for regions of the network associated with the phenotype using [SConES](https://academic.oup.com/bioinformatics/article/29/13/i171/198210) or [SigMod](https://academic.oup.com/bioinformatics/article/33/10/1536/2874362).

# Installation

Install `martini` like any other Bioconductor package:

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

# 3. Run SConES, finding the best parameters by cross-validation
res <- scones.cv(minigwas, gs)

# the output is an igraph subnetwork containing the selected SNPs
res
# IGRAPH d9128a0 UNW- 12 10 -- 
# + attr: name (v/c), chr (v/n), pos (v/n), weight (e/n)
# + edges from d9128a0 (vertex names):
# [1] 1A1--1A2 1A2--1A3 1A3--1A4 1A4--1A5 1A5--1A6 2C1--2C2 2C2--2C3 2C3--2C4 2C4--2C5 2C5--2C6
```

Please, refer to the vignettes for more detailed usage examples. `martini` results can be further examined using the [blur](https://github.com/hclimente/blur) package.

# Citation

A more detailed description can be found in [the pre-print](https://www.biorxiv.org/content/10.1101/2021.01.25.428047v1). If you use `martini` in your work, please cite us:

```
@article{martini2021,
	title = {martini: an {R} package for genome-wide association studies using {SNP} networks},
  author = {Climente-González, Héctor and Azencott, Chloé-Agathe},
	url = {http://biorxiv.org/lookup/doi/10.1101/2021.01.25.428047},
	journal = {bioRxiv},
	month = jan,
	year = {2021},
	doi = {10.1101/2021.01.25.428047}
}
```
