#' Measure epistasis between SNPs
#' 
#' @description Measure epistasis between SNPs along the edges of net.
#' 
#' @param gwas A GWAS experiment.
#' @param net A network of SNPs.
#' @return A dataframe with the chi-squared score for each pair of tested SNPs.
#' @importFrom igraph as_data_frame
#' @importFrom stats chisq.test
#' @export
measure_epistasis <- function(gwas, net) {
  
  edges <- as_data_frame(net, what = "edges")
  colnames(edges) <- c("snp1", "snp2")
  
  X <- as(gwas$genotypes, "numeric")
  colnames(X) <- as.character(gwas$map[,2])
  
  edges$score <- mapply(function(x, y) {
    chisq.test(X[,x], X[,y])$statistic
  }, edges$snp1, edges$snp2)
  
  return(edges)
  
}

#' Incorporate epistasis into the network
#' 
#' @description Incorporate epistatic interactions in the network.
#' 
#' @param net A network of SNPs.
#' @param scores A dataframe with pairs of SNPs and the measure of the epistatic effect.
#' @return A network where the weight of the edges reflect the strength of the epistatic interaction.
#' @importFrom igraph E set_edge_attr %>%
#' @export
weight_epistasis_network <- function(net, scores) {
  
  colnames(scores) <- c("snp1", "snp2", "score")
  snps <- subset(scores, select = c("snp1", "snp2")) %>% t %>% c
  
  net <- set_edge_attr(net, "weight", index = E(net, P=snps), value = scores$score)
  
  return(net)
  
}