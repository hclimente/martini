#' Get GS network.
#' 
#' @description Creates a network of SNPs where each SNP is connected to its adjacent SNPs in the genome sequence. Corresponds to the GS network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @return An igraph network of the GS network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171-179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @export
get_GS_network <- function(gwas)  {
  
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  map <- subset(map, select = c("chr","pos","snp"))
  map <- unique(map)
  map <- map[with(map, order(chr,pos)),]
  
  gs <- by(map, map$chr, function(x){
    chr <- unique(x$chr)
    snp1 <- head(x$snp, n = length(x$snp) - 1)
    snp2 <- tail(x$snp, n = length(x$snp) - 1)
    data.frame(snp1 = snp1, snp2 = snp2)
  })
  gs <- do.call("rbind", gs)
  gs <- graph_from_data_frame(gs, directed = FALSE)
  gs <- simplify(gs)
  
  gs <- set_vertex_attr(gs, "chr", index = match(map$snp, V(gs)$name), map$chr)
  gs <- set_vertex_attr(gs, "pos", index = match(map$snp, V(gs)$name), map$pos)

  return(gs)
  
}