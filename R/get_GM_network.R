#' Get GM network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the GS network and, in addition, to all the other SNPs pertaining to the same gene. Corresponds to the GM network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snp2gene A data frame with two columns: snp id (1st column) and gene it maps to (2nd column).
#' @return An igraph network of the GM network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171–179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @export
get_GM_network <- function(gwas, snp2gene)  {
  
  colnames(snp2gene) <- c("snp","gene")
  
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  map <- merge(map, snp2gene)
  map <- subset(map, select = c("snp","gene"))
  # in some cases the same position is linked to two different variants
  map <- unique(map)
  
  gm <- by(map, map$gene, function(x){
    if (nrow(x) > 1){
      comb <- combn(x$snp, 2)
      data.frame(snp1 = comb[1,], snp2 = comb[2,])
    }
  })
  gm <- do.call("rbind", gm)
  gm <- graph_from_data_frame(gm, directed = FALSE)
  gs <- get_GS_network(gwas)
  gm <- simplify(gm + gs)

  return(gm)
}