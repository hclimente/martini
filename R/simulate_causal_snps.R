#' Simulate causal SNPs
#' 
#' @description Takes the biggest clique of the network, and randomly picks n SNPs from it.
#' 
#' @param net An igraph network that connects the SNPs.
#' @param n Number of causal SNPs to return.
#' @return A boolean vector with as many elements as SNPs.
#' @export
simulate_causal_snps <- function(gwas, net, n) {

  idx <- 1
  
  cliques <- largest_cliques(net)
  myClique <- cliques[[idx]]
  
  # randomly select snps
  causalIds <- sample(myClique, n)
  
  if ("snp.names" %in% colnames(gwas$map))
    causal <- gwas$map$snp.names %in% names(causalIds)
  else if ("snp" %in% colnames(gwas$map))
    causal <- gwas$map$snp %in% names(causalIds)
  
  return(causal)
}