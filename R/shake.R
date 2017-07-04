#' GWAS incorporating networks.
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being connected in an underlying network (Azencott et al., 2013).
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param net An igraph network that connects the SNPs.
#' @param ... Extra arguments for shake.
#' @return A copy of the SnpMatrix$map object, with the following additions:
#' \itemize{
#' \item{C: contains the univariate association score for every single SNP.}
#' \item{selected: logical vector indicating if the SNP was selected by shake or not.}
#' }
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171–179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @export
shake <- function(gwas, net, ...) {

  X <- as(gwas$genotypes, "numeric")
  Y <- gwas$fam$affected
  
  # remove redundant edges and self-edges
  net <- simplify(net)
  W <- get.adjacency(net, type="both", sparse = TRUE)
  
  # order according to order in map
  W <- W[gwas$map$snp.names, gwas$map$snp.names]
  
  settings <- get_shake_settings(...)
  
  gin <- run_shake(X, Y, W, settings)
  cat("eta =", gin$eta, "\nlambda =", gin$lambda)
  
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  map$C <- gin$scores
  map$selected <- as.logical(gin$indicator)
  
  return(map)
  
}