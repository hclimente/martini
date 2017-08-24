#' Calculate an empirical p-value for \code{find_cones} results.
#' 
#' @description Do a permutation-based test to assess the statistical significance of each of the clusters obtained through shake. For a cluster of size k, k interconnected SNPs are picked N times, and their joint association score is calculated to come up with an estimation of the distribution.
#' 
#' @param cones Results from \code{find_cones}.
#' @param net The same SNP network provided to \code{find_cones}.
#' @param N Integer with the name of permutations.
#' @return An empirical p-value for each of the SNP clusters. Please, note that the minimum possible p-value from an empirical distribution is set to 1/(N+1). Clusters composed by a single SNP will not have a empirical p-value.
#' @importFrom stats ecdf
#' @importFrom igraph vcount random_walk
#' @export
test_cones_clusters <- function(cones, net, N = 100000) {
  
  numSNPs <- vcount(net)
  selected <- subset(cones, selected)
  
  # calculate one ecdf for each cluster size
  clusterSizes <- unique(table(selected$cluster))
  names(clusterSizes) <- clusterSizes
  
  ecdfs <- lapply(clusterSizes, function(n) {
    if (n > 1){
      sampled_C <- lapply(1:N, function(i){
        v <- sample(1:numSNPs, 1)
        snpCluster_i <- random_walk(net, v, n)
        cones_cluster <- subset(cones, snp %in% snpCluster_i)
        sum(cones_cluster$c)
      })
      sampled_C <- do.call("c", sampled_C)
      
      ecdf(sampled_C)
    }
  })
  
  clusters <- by(selected, selected$cluster, function(k){
    n <- nrow(k)
    i <- unique(k$cluster)
    
    if (n > 1){
      snpCluster_C <- sum(k$c)
      p <- 1 - ecdfs[[as.character(n)]](snpCluster_C)
    } else {
      p <- NA
    }
    
    data.frame(k = i, ncomponents = n, p = p)
  })
  
  clusters <- do.call("rbind", clusters)
  # minimum p-vaue = 1 / (N + 1)
  clusters$p[clusters$p == 0] = 1 / (N + 1)
  
  return(clusters)

}
