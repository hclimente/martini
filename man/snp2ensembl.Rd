% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db.R
\name{snp2ensembl}
\alias{snp2ensembl}
\title{Map SNPs to Ensembl genes.}
\usage{
snp2ensembl(gwas, organism = 9606, flank = 0)
}
\arguments{
\item{gwas}{A SnpMatrix object with the GWAS information.}

\item{organism}{Tax ID of the studied organism. The default is 9606 (human).}

\item{flank}{A number with the flanking regions around genes to be considered
part of the gene i.e. SNPs mapped to them will be considered mapped to the
gene.}
}
\value{
A data.frame with two columns: one for the SNP and another for the
gene it has been mapped to.
}
\description{
Maps SNPs from a GWAS experiment to genes.
}
\keyword{internal}
