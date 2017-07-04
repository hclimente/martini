#' Simulate causal SNPs
#' 
#' @description Takes the biggest clique of the network, and randomly picks n SNPs from it.
#' 
#' @param net An igraph network that connects the SNPs.
#' @param n Number of causal SNPs to return.
#' @return A boolean vector with as many elements as SNPs.
#' @export
simulate_causal_snps <- function(net, n) {

  idx <- 1
  
  graph <- graph_from_adjacency_matrix(net, mode = "undirected", diag = FALSE)
  cliques <- largest_cliques(graph)
  cliqueIdx <- cliques[[idx]]
  
  # randomly select snps
  causalSnpIdx <- sample(as.numeric(cliqueIdx), n)
  
  causalSnp <- rep(FALSE, nrow(net))
  causalSnp[causalSnpIdx] <- TRUE
  
  return(causalSnp)
}