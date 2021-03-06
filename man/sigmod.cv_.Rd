% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nongwas.R
\name{sigmod.cv_}
\alias{sigmod.cv_}
\title{Find connected explanatory features}
\usage{
sigmod.cv_(X, y, featnames, net)
}
\arguments{
\item{X}{n x d design matrix}

\item{y}{Vector of length n with the outcomes}

\item{featnames}{Vector of length d with the feature names}

\item{net}{An igraph network that connects the SNPs.}
}
\value{
A copy of the \code{SnpMatrix$map} \code{data.frame}, with the 
following additions:
\itemize{
\item{c: contains the univariate association score for every single SNP.}
\item{selected: logical vector indicating if the SNP was selected by SConES 
or not.}
\item{module: integer with the number of the module the SNP belongs to.}
}
}
\description{
Finds the features maximally associated with a phenotype while 
being connected in an underlying network. Select the hyperparameters by
cross-validation.
}
\examples{
X <- as(minigwas[['genotypes']], 'numeric')
X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))
gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
sigmod.cv_(X, minigwas[['fam']]$affected, minigwas[['map']]$snp, gi)
}
