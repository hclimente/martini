#' Include LD information in the network.
#' 
#' @description Include LD information in the SNP network. Only work for human
#' SNPs.
#' 
#' @param net A SNP network.
#' @param ld Data frame with linkage disequilibrium, like \code{snpStats::ld} 
#' output.
#' @param method How to incorporate LD values into the network.
#' @return An SNP network where the edges weight 1 - LD, measured as Pearson
#' correlation.
#' @importFrom igraph E %>% set_edge_attr delete_edges get.edgelist
#' @importFrom Matrix summary
#' @importFrom stats cor
#' @examples 
#' ld <- snpStats::ld(minigwas$genotypes, depth = 2, stats = "R.squared")
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' ldGi <- ldweight_edges(gi, ld)
#' @export
ldweight_edges <- function(net, ld, method = "inverse") {
  
  edges <- as.data.frame(get.edgelist(net))
  edges <- paste(edges[,1], edges[,2], sep = "-")
  ldDf <- Matrix::summary(ld)
  ldDf$key <- paste(rownames(ld)[ldDf$i], rownames(ld)[ldDf$j], sep = "-")
  ldDf <- subset(ldDf, key %in% edges)
  idx <- match(ldDf$key, edges)

  if (method == "inverse") {
    net <- set_edge_attr(net, "weight", index = idx, value = 1 / (1 + ldDf$x))
  } else if (method == "subtraction") {
    net <- set_edge_attr(net, "weight", index = idx, value = 1 - ldDf$x)
  }
  
  if (any(is.na(E(net)$weight))) {
    stop("NA values as edge weights.")
  } else if (any(E(net)$weight < 0)) {
    stop("Edge weights cannpt be negative.")
  }
  
  # remove edges with 0 weight
  zeroE <- E(net)[E(net)$weight == 0]
  net <- delete_edges(net, zeroE)
  
  return(net)
  
}

#' LD-remove redundant SNPs
#' 
#' @description Removes SNPs belonging to the same LD-block, leaving only one. 
#' An LD-block is defined as sets of SNPs that have a LD measure with the 
#' adjacent SNP higher than the cut-off.
#' 
#' @param gwas A GWAS experiment.
#' @param ld A matrix of # SNPs x # SNPs representing the LD between pairs of 
#' SNPs. For this function to work, only the LD with the adjacent SNPs is 
#' required.
#' @param cutoff Minimum LD to be considered part of the same LD-block.
#' @return A copy of the GWAS experiment without redundant SNPs.
#' @examples 
#' ld <- snpStats::ld(minigwas$genotypes, depth = 2, stats = "R.squared")
#' ld_prune(minigwas, ld)
#' @export
ld_prune <- function(gwas, ld, cutoff = 0.8) {
  
  if (is_coherent(gwas) & any(gwas$map[,2] != colnames(ld))) {
    stop("gwas$map and ld SNP order differ.")
  }
  
  map <- gwas$map
  
  b <- 1
  blocks <- numeric(nrow(ld))
  blocks[1] <- 1
  
  for (i in 2:nrow(ld)) {
    thisLd <- ld[i - 1,i]
    if(is.na(thisLd) | thisLd < cutoff) {
      b <- b + 1
    }
    blocks[i] <- b
    
  }

  representative <- by(map, blocks, function(block) {
    
    snps <- as.character(block[,2])
    id <- ceiling(length(snps)/2)
    
    return(snps[id])
    
  }) %>% as.character
  
  gwas$map <- map[map[,2] %in% representative, ]
  gwas$genotypes <- gwas$genotypes[, representative]
  gwas$pruning <- data.frame(snp = map[,2], block = blocks)
  
  return(gwas)
  
}