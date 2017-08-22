#' Return groups of interconnected SNPs.
#' 
#' @description Find clusters composed by interconnected SNPs.
#' 
#' @param map Results from \code{shake}.
#' @param net The same SNP network provided to \code{shake}.
#' @return A list with the clusters of selected SNPs.
#' @export
cluster_snps <- function(map, net) {
  
  selected <- subset(map, selected)
  subnet <- igraph::induced_subgraph(net, selected$snp)
  
  clusters <- components(subnet)
  map$cluster <- NA
  map$cluster[map$selected] <- clusters$membership[order(match(names(clusters$membership), map$snp))]
  
  return(map)
}