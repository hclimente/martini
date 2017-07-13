#' Simulate causal SNPs
#' 
#' @description Selects randomly a set of loosely interconnected SNPs
#' 
#' @param net An igraph network that connects the SNPs.
#' @param n Number of causal SNPs to return.
#' @return A boolean vector with as many elements as SNPs.
#' @export
simulate_causal_snps <- function(gwas, net, n) {

  k <- estimate_closeness(net, sample(V(net), floor(length(V(net)) * 0.05)), cutoff=20)
  
  seed <- names(tail(sort(k), n=1))
  x <- lapply(seq(100), function(x) {
    names(random_walk(net, seed, 2*n) )
  })
  
  causalIds <- names(head(sort(table(do.call("c", x)), decreasing=T), n=n))

  return(causalIds)
}