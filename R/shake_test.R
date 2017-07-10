#' Calculate an empirical p-value for shake results.
#' 
#' @description Do a permutation-based test to assess the statistical significance of each of the clusters obtained through shake.
#' 
#' @param map Results from \code{shake}.
#' @param net The same SNP network provided to \code{shake}.
#' @param n Integer with the name of permutations.
#' @return An empirical p-value for each of the SNP clusters.
#' @export
shake_test <- function(map, net, n = 100000) {
  
  numSNPs <- vcount(net)
  selected <- subset(map, selected)
  
  clusters <- by(map, map$cluster, function(k){
    
    snpCluster_C <- sum(k$C)
    
    sampled_C <- lapply(1:n, function(i){
      v <- sample(1:numSNPs, 1)
      snpCluster_i <- random_walk(net, v, nrow(k))
      map_cluster <- subset(map, snp %in% snpCluster_i)
      sum(map_cluster$C)
    })
    sampled_C <- do.call("c", sampled_C)
    
    C_df <- ecdf(sampled_C)
    data.frame(k = unique(k$cluster), ncomponents = nrow(k), p = C_df(snpCluster_C))
  })
  
  clusters <- do.call("rbind", clusters)
  
  return(clusters)

}