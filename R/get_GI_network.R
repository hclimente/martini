#' Get GI network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the GM network and, in addition, to all the other SNPs pertaining to any interactor of the gene it is mapped to. Corresponds to the GI network described by Azencott et al.
#' 
#' @param map A map object with the SNP position information.
#' @param snp2gene A data frame with two columns: snp and gene. The first column contains the id of the SNP; the second, the gene which is mapped to.
#' @param tab A data frame describing protein-protein interactions in TAB format.
#' @return An igraph network of the GI network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171–179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @export
get_GI_network <- function(map, snp2gene, tab)  {
  
  genes <- read.delim2(snp2gene)
  colnames(genes) <- c("snp","gene")
  
  map <- read.delim2(map, sep = " ", header = FALSE)
  colnames(map) <- c("chr","snp","cm","pos")
  map <-  subset(map, select = c("chr","snp","pos"))
    
  ppi <- read.delim2(tab)
  ppi$gene1 <- ppi$OFFICIAL_SYMBOL_FOR_A
  ppi$gene2 <- ppi$OFFICIAL_SYMBOL_FOR_B
  ppi <- subset(ppi, gene1 != gene2, select = c("gene1","gene2"))
  ppi <- unique(ppi)
  
  gi <- merge(ppi, genes, by.x = "gene1", by.y = "gene")
  gi <- merge(gi, genes, by.x = "gene2", by.y = "gene")
  gi <- merge(gi, map, by.x = "snp.x", by.y = "snp")
  gi <- merge(gi, map, by.x = "snp.y", by.y = "snp")
  gi$chr1 <- gi$chr.x
  gi$pos1 <- gi$pos.x
  gi$chr2 <- gi$chr.y
  gi$pos2 <- gi$pos.y
  gi <- subset(gi, chr1 != chr2 | pos1 != pos2, select = c("chr1", "pos1", "chr2", "pos2"))
  
  gi <- graph_from_data_frame(gi, directed = FALSE)
  gi_adj <- get.adjacency(gi, type="both", sparse = TRUE)
  
  gm_adj <- getGMNet(map, snp2gene)
  
  return(gi_adj + gm_adj)
  
}