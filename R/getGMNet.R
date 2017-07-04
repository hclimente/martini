getGMNet <- function(map, snp2gene)  {
  
  genes <- read.delim2(snp2gene)
  colnames(genes) <- c("snp","gene")
  
  map <- read.delim2(map, sep = " ", header = FALSE)
  colnames(map) <- c("chr","snp","cm","pos")
  map <- merge(map, genes)
  map <- subset(map, select("chr","gene","pos"))
  # in some cases the same position is linked to two different variants
  map <- unique(map)
  
  gm <- by(map, map$gene, function(x){
    chr <- unique(x$chr)
    if (nrow(x) > 1){
      comb <- combn(x$pos, 2)
      data.frame(chr1 = chr, pos1 = comb[1,], chr2 = chr, pos2 = comb[2,])
    }
  })
  gm <- do.call("rbind", gm)
  gm <- graph_from_data_frame(gm, directed = FALSE)
  gm_adj <- get.adjacency(gm, type="both", sparse = TRUE)
  
  gs_adj <- getGSNet(map)

  return(gm_adj + gs_adj)
}