get_GS_network <- function(map)  {
  
  map <- read.delim2(map, sep = " ", header = FALSE)
  colnames(map) <- c("chr","gene","cm","pos")
  map <- subset(map, select = c("chr","pos"))
  map <- unique(map)
  map <- map[with(map, order(chr,pos)),]
  
  gs <- by(map, map$chr, function(x){
    chr <- unique(x$chr)
    pos1 <- head(x$pos, n = length(x$pos) - 1)
    pos2 <- tail(x$pos, n = length(x$pos) - 1)
    data.frame(snp1 = paste(chr, pos1, sep = "_"), 
               snp2 = paste(chr, pos2, sep = "_"))
  })
  gs <- do.call("rbind", gs)
  gs <- graph_from_data_frame(gs, directed = FALSE)
  gs_adj <- get.adjacency(gs, type="both", sparse = TRUE)
  
  return(gs_adj)
  
}