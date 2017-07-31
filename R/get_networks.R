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
  gs <- igraph::graph_from_data_frame(gs, directed = FALSE)
  gs <- igraph::simplify(gs)
  
  gs <- igraph::set_vertex_attr(gs, "chr", index = match(map$snp, V(gs)$name), map$chr)
  gs <- igraph::set_vertex_attr(gs, "pos", index = match(map$snp, V(gs)$name), map$pos)

  return(gs)
  
}

#' Get GM network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the GS network and, in addition, to all the other SNPs pertaining to the same gene. Corresponds to the GM network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snp2gene A data frame with two columns: snp id (1st column) and gene it maps to (2nd column).
#' @return An igraph network of the GM network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171-179. \url{https://doi.org/10.1093/bioinformatics/btt238}
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
  
  gs <- get_GS_network(gwas)
  
  if (! is.null(gm)) {
    gm <- igraph::graph_from_data_frame(gm, directed = FALSE)
    gm <- igraph::simplify(gm + gs)
    gm <- igraph::set_vertex_attr(gm, "gene", index = match(map$snp, V(gm)$name), map$gene)
  } else {
    warning("insufficient information to add gene information")
    gm <- gs
  }
  
  return(gm)
}

#' Get GI network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the GM network and, in addition, to all the other SNPs pertaining to any interactor of the gene it is mapped to. Corresponds to the GI network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snp2gene A data frame with minimum two columns: snp id (1st column) and gene it maps to (2nd column).
#' @param ppi A data frame describing protein-protein interactions with at least two colums. The first two columns must be the gene ids of the interacting proteins.
#' @return An igraph network of the GI network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus association mapping with graph cuts. Bioinformatics, 29(13), 171-179. \url{https://doi.org/10.1093/bioinformatics/btt238}
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
  
  if (nrow(snp2snp) == 0)
    warning("no matches between genes in snp2gene and PPI. No information about PPI will be added.")
  
  gi <- igraph::graph_from_data_frame(snp2snp, directed = FALSE)
  gm <- get_GM_network(gwas, snp2gene)
  gi <- igraph::simplify(gm + gi)
  
  return(gi)
  
}