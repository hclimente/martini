#' Simulate causal SNPs
#' 
#' @description Selects randomly a clique of interconnected SNPs. If the SNP network contains a "gene" vertex attribute, it tries to pick SNPs from the same gene and from, at least, one interactor. Else, it picks SNPs from the largest clique in the network.
#' 
#' @param net An igraph network that connects the SNPs.
#' @param n Number of causal SNPs to return.
#' @return A vector with the ids of the simulated causal SNPs.
#' @export
simulate_causal_snps <- function(net, n) {
  
  if (! is.null(igraph::vertex_attr(net, "gene"))) {
    genes <- names(which(table(igraph::V(net)$gene) > 1))
    
    repeat {
      g <- sample(genes, 1)
      seed <- igraph::V(net)[which(igraph::V(net)$gene == g)][1]
      
      neighboringGenes <- igraph::neighbors(net, seed)$gene %>% na.omit %>% unique
      neighboringGenes <- intersect(genes, neighboringGenes)
      neighbors <- igraph::V(net)$gene %in% neighboringGenes %>% which %>% igraph::V(net)[.]
      
      if ( length(neighbors) >= n & any(neighboringGenes != g) ) {
        causal <- sample(neighbors, n)
        genesInvolved <- unique(causal$gene)
        
        if (length(genesInvolved) > 1)
          break
      }
    }
    
  } else {
    largestClique <- igraph::largest_cliques(net)[[1]]
    causal <- sample(largestClique, n)
  }
  
  return(causal)
}