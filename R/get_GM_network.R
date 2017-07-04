#' Get GM network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the GS network and, in addition, to all the other SNPs pertaining to the same gene. Corresponds to the GM network described by Azencott et al.
#' 
#' @param map A map object with the SNP position information.
#' @param snp2gene A data frame with two columns: snp and gene. The first column contains the id of the SNP; the second, the gene which is mapped to.
#' @return An igraph network of the GM network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171–179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @export
get_GM_network <- function(map, snp2gene)  {
  
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