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