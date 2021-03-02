#' Get genomic sequence network
#' 
#' @description Creates a network of SNPs where each SNP is connected to its 
#' adjacent SNPs in the genome sequence. Corresponds to the genomic sequence 
#' (GS) network described by Azencott et al.
#' @template params_gwas
#' @return An igraph network of the GS network of the SNPs.
#' @template reference_azencott
#' @importFrom igraph graph_from_data_frame simplify set_vertex_attr V
#' set_edge_attr
#' @importFrom stats aggregate
#' @importFrom utils combn head tail
#' @examples
#' get_GS_network(minigwas)
#' @export
get_GS_network <- function(gwas)  {
  
  map <- sanitize_map(gwas)
  map <- subset(map, select = c("chr","pos","snp"))
  map <- unique(map)
  map <- map[with(map, order(chr,pos)),]
  
  gs <- by(map, map[,'chr'], function(x){
    chr <- unique(x[,'chr'])
    snp1 <- head(x[,'snp'], n = length(x[,'snp']) - 1)
    snp2 <- tail(x[,'snp'], n = length(x[,'snp']) - 1)
    data.frame(snp1 = snp1, snp2 = snp2)
  })
  gs <- do.call(rbind, gs)
  gs <- graph_from_data_frame(gs, directed = FALSE)
  gs <- simplify(gs)
  
  gs <- set_vertex_attr(gs, "chr", 
                        index = match(map[,'snp'], V(gs)$name), map[,'chr'])
  gs <- set_vertex_attr(gs, "pos", 
                        index = match(map[,'snp'], V(gs)$name), map[,'pos'])
  gs <- set_edge_attr(gs, "weight", value = 1)

  return(gs)
  
}

#' Get gene membership network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the
#' \link[=get_GS_network]{GS} network and, in addition, to all the other SNPs 
#' pertaining to the same gene. Corresponds to the gene membership (GM) network 
#' described by Azencott et al.
#' @template params_gwas
#' @template params_organism
#' @template params_snpMapping
#' @template params_col_genes
#' @return An igraph network of the GM network of the SNPs.
#' @template reference_azencott
#' @importFrom igraph %>% graph_from_data_frame simplify set_vertex_attr V
#' set_edge_attr
#' @importFrom utils combn
#' @examples 
#' get_GM_network(minigwas, snpMapping = minisnpMapping)
#' @export
get_GM_network <- function(gwas, organism = 9606, 
                           snpMapping = snp2ensembl(gwas, organism),
                           col_genes = c('snp','gene')) {
  
  snpMapping <- sanitize_snpMapping(snpMapping, col_genes) 
  map <- sanitize_map(gwas)
  
  gs <- get_GS_network(gwas)
  
  map <- merge(map, snpMapping)
  map <- subset(map, select = c("snp","gene"))
  # in some cases the same position is linked to two different variants
  map <- unique(map)
  map[,'snp'] <- as.character(map[,'snp'])
  map[,'gene'] <- as.character(map[,'gene'])
 
  # use only genes that map to more than 1 snp
  gene_freq <- table(map[,'gene'])
  genes <- names(gene_freq)[gene_freq > 1]
   
  if (length(genes)) {
    
    map <- subset(map, gene %in% genes)
    
    gm <- do.call(cbind, tapply(map[,'snp'], map[,'gene'], combn, 2)) %>% 
      t %>% 
      as.data.frame
    gm <- graph_from_data_frame(gm, directed = FALSE)
    gm <- simplify(gm + gs)
    gm <- set_vertex_attr(gm, 'gene', index = match(map[,'snp'], V(gm)$name), 
                          map[,'gene'])
    
    nGenes <- aggregate(gene ~ ., data=map, length)
    gm <- set_vertex_attr(gm, "nGenes", value = 0)
    gm <- set_vertex_attr(gm, "nGenes", 
                          index = match(nGenes[,'snp'], V(gm)$name),nGenes$gene)
  } else {
    warning("insufficient information to add gene information")
    gm <- gs
  }
  
  gm <- set_edge_attr(gm, "weight", value = 1)
  
  return(gm)
}

#' Get gene-interaction network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the
#' \link[=get_GM_network]{GM} network and, in addition, to all the other SNPs 
#' pertaining to any interactor of the gene it is mapped to. Corresponds to the 
#' gene-interaction (GI) network described by Azencott et al.
#' 
#' @template params_gwas
#' @template params_organism
#' @template params_snpMapping
#' @param ppi A data.frame describing protein-protein interactions with at least
#' two colums. Gene ids must be the contained in snpMapping. Unless column names
#' are specified using \code{col_ppi}, involved columns must be named 
#' \code{gene1} and \code{gene2}.
#' @template params_col_genes
#' @param col_ppi Optional, length-2 character vector with the names of the two 
#' columns involving the protein-protein interactions.
#' @template params_flush
#' @return An igraph network of the GI network of the SNPs.
#' @template reference_azencott
#' @importFrom igraph graph_from_data_frame simplify set_edge_attr
#' @examples 
#' get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' @export
get_GI_network <- function(gwas, organism = 9606,
                           snpMapping = snp2ensembl(gwas, organism),
                           ppi = get_gxg('biogrid', organism, flush),
                           col_ppi = c('gene1','gene2'),
                           col_genes = c('snp','gene'),
                           flush = FALSE) {
  
  snpMapping <- sanitize_snpMapping(snpMapping, col_genes) 
  
  map <- sanitize_map(gwas)
  map <-  subset(map, select = c("chr","snp","pos"))
  
  # read tab file  
  ppi <- subset(ppi, select = col_ppi)
  colnames(ppi) <- c("gene1", "gene2")
  ppi <- unique(ppi)
  # remove self-interactions
  ppi <- subset(ppi, gene1 != gene2)
  
  # match all SNPs to pairwise PPI
  snpMapping <- subset(snpMapping, select = col_genes)
  colnames(snpMapping) <- c('snp','gene')
  snp2snp <- merge(ppi, snpMapping, by.x = "gene1", by.y = "gene")
  snp2snp <- merge(snp2snp, snpMapping, by.x = "gene2", by.y = "gene")
  snp2snp <- merge(snp2snp, map, by.x = "snp.x", by.y = "snp")
  snp2snp <- merge(snp2snp, map, by.x = "snp.y", by.y = "snp")
  
  if (nrow(snp2snp) == 0) {
    warning("no matches between genes in snpMapping and PPI. No information about PPI will be added.")
  }
  
  gi <- graph_from_data_frame(snp2snp, directed = FALSE)
  gm <- get_GM_network(gwas, snpMapping=snpMapping)
  gi <- simplify(gm + gi)
  
  gi <- set_edge_attr(gi, "weight", value = 1)
  
  return(gi)
  
}
