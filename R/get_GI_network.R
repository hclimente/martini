get_GI_network <- function(tab, snp2gene, map)  {
  
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