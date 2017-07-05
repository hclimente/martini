#' Get GI network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the GM network and, in addition, to all the other SNPs pertaining to any interactor of the gene it is mapped to. Corresponds to the GI network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snp2gene A data frame with minimum two columns: snp id (1st column) and gene it maps to (2nd column).
#' @param ppi A data frame describing protein-protein interactions with at least two colums. The first two columns must be the gene ids of the interacting proteins.
#' @return An igraph network of the GI network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171–179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @export
get_GI_network <- function(gwas, snp2gene, ppi)  {
  
  colnames(snp2gene) <- c("snp","gene")
  
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  map <-  subset(map, select = c("chr","snp","pos"))
  
  # read tab file  
  ppi <- ppi[,c(1,2)]
  colnames(ppi) <- c("gene1", "gene2")
  ppi <- unique(ppi)
  # remove self-interactions
  ppi <- subset(ppi, gene1 != gene2)
  
  # match all SNPs to pairwise PPI
  snp2snp <- merge(ppi, snp2gene, by.x = "gene1", by.y = "gene")
  snp2snp <- merge(snp2snp, snp2gene, by.x = "gene2", by.y = "gene")
  snp2snp <- merge(snp2snp, map, by.x = "snp.x", by.y = "snp")
  snp2snp <- merge(snp2snp, map, by.x = "snp.y", by.y = "snp")

  # remove self-interactions
  snp2snp <- subset(snp2snp, chr.x != chr.y & pos.x != pos.y, select = c("snp.x", "snp.y"))
  
  gi <- graph_from_data_frame(snp2snp, directed = FALSE)
  gm <- get_GM_network(gwas, snp2gene)
  gi <- simplify(gm + gi)
  
  return(gi)
  
}