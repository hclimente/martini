#' Get genomic sequence network.
#' 
#' @description Creates a network of SNPs where each SNP is connected to its 
#' adjacent SNPs in the genome sequence. Corresponds to the genomic sequence 
#' (GS) network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @return An igraph network of the GS network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., &
#' Borgwardt, K. M. (2013). Efficient network-guided multi-locus association
#' mapping with graph cuts. Bioinformatics, 29(13), 171-179.
#' \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph graph_from_data_frame simplify set_vertex_attr V
#' set_edge_attr
#' @importFrom stats aggregate
#' @importFrom utils combn head tail
#' @examples
#' get_GS_network(minigwas)
#' @export
get_GS_network <- function(gwas)  {
  
  map <- gwas[["map"]]
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
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
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Tax ID of the studied organism. Required if snpMapping is not
#' provided.
#' @param snpMapping A data.frame informing how SNPs map to genes. It contains 
#' minimum two columns: SNP id and a gene it maps to. Each row corresponds to 
#' one gene-SNP mapping. Unless column names are specified using 
#' \code{col_genes}, involved columns must be named \code{'snp'} and 
#' \code{'gene'}.
#' @param col_genes Optional, length-2 character vector with the names of the 
#' two columns involving the SNP-gene mapping. The first element is the column 
#' of the SNP, and the second is the column of the gene.
#' @return An igraph network of the GM network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., &
#' Borgwardt, K. M. (2013). Efficient network-guided multi-locus association
#' mapping with graph cuts. Bioinformatics, 29(13), 171-179.
#' \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph delete_edge_attr delete_vertex_attr simplify union
#' @examples 
#' get_GM_network(minigwas, snpMapping = minisnpMapping)
#' @export
get_GM_network <- function(gwas, organism = 9606, 
                           snpMapping = snp2gene(gwas, organism),
                           col_genes = c('snp','gene')) {
    
    gs <- get_GS_network(gwas)
    ft <- get_feature_network(gwas, snpMapping, col_genes, 'gene')
    
    gm <- union(ft, gs)
    
    gm <- set_vertex_attr(gm, "nGenes", 
                          value = ifelse(is.na(V(gm)$n_features), 0, V(gm)$n_features))
    gm <- delete_vertex_attr(gm , 'n_features')
    
    gm <- set_edge_attr(gm, "weight", value = 
                        mapply(max, E(gm)$weight_1, E(gm)$weight_2, na.rm = T))
    gm <- delete_edge_attr(gm , 'weight_1')
    gm <- delete_edge_attr(gm , 'weight_2')
    
    return(gm)

}

#' Get gene-interaction network.
#' 
#' @description Creates a network of SNPs where each SNP is connected as in the
#' \link[=get_GM_network]{GM} network and, in addition, to all the other SNPs 
#' pertaining to any interactor of the gene it is mapped to. Corresponds to the 
#' gene-interaction (GI) network described by Azencott et al.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Tax ID of the studied organism. Required if snpMapping is not
#' provided.
#' @param snpMapping A data.frame informing how SNPs map to genes. It contains 
#' minimum two columns: SNP id and a gene it maps to. Each row corresponds to 
#' one gene-SNP mapping. Unless column names are specified using 
#' \code{col_genes}, involved columns must be named \code{'snp'} and 
#' \code{'gene'}.
#' @param ppi A data.frame describing protein-protein interactions with at least
#' two colums. Gene ids must be the contained in snpMapping. Unless column names
#' are specified using \code{col_ppi}, involved columns must be named 
#' \code{gene1} and \code{gene2}.
#' @param col_genes Optional, length-2 character vector with the names of the 
#' two columns involving the SNP-gene mapping. The first element is the column 
#' of the SNP, and the second is the column of the gene.
#' @param col_ppi Optional, length-2 character vector with the names of the two 
#' columns involving the protein-protein interactions.
#' @return An igraph network of the GI network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., &
#' Borgwardt, K. M. (2013). Efficient network-guided multi-locus association
#' mapping with graph cuts. Bioinformatics, 29(13), 171-179.
#' \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph graph_from_data_frame simplify set_edge_attr
#' @examples 
#' get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' @export
get_GI_network <- function(gwas, organism, 
                           snpMapping = snp2gene(gwas, organism), 
                           ppi = get_ppi(organism), 
                           col_ppi = c('gene1','gene2'),
                           col_genes = c('snp','gene')) {
  
  colnames(snpMapping) <- c("snp","gene")
  
  map <- gwas[["map"]]
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
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

#' Get feature membership network.
#' 
#' @description Creates a network of SNPs where each SNP is connected to all the
#' other SNPs mapped to the same feature.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snpMapping A data.frame informing how SNPs map to features It contains
#' at least two columns: SNP id and a feature it maps to. Each row corresponds 
#' to one gene-SNP mapping. Unless column names are specified using 
#' \code{col_features}, involved columns must be named \code{'snp'} and 
#' \code{'feature'}.
#' @param col_features Optional, length-2 character vector with the names of the 
#' two columns involving the SNP-gene mapping. The first element is the column 
#' of the SNP, and the second is the column of the gene.
#' @param feature_type Feature type.
#' @return An igraph network where SNPs mapping to the same feature form 
#' cliques.
#' @importFrom igraph %>% graph_from_data_frame simplify set_vertex_attr V
#' set_edge_attr
#' @importFrom utils combn
#' @export
get_feature_network <- function(gwas, snpMapping,
                                col_features = c('snp','feature'),
                                feature_type = 'feature') {
    
    snpMapping <- subset(snpMapping, select = col_features)
    colnames(snpMapping) <- c('snp','feature')
    
    map <- gwas[["map"]]
    colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
    map <- merge(map, snpMapping)
    map <- subset(map, select = c("snp","feature"))
    # in some cases the same position is linked to two different variants
    map <- unique(map)
    map[,'snp'] <- as.character(map[,'snp'])
    map[,'feature'] <- as.character(map[,'feature'])
    
    nSnps <- aggregate(snp ~ ., data=map, length)
    map <- subset(map, feature %in% nSnps[nSnps$snp > 1, 'feature'])
    
    if (nrow(map) > 0) {
        ft <- do.call(cbind, tapply(map[,'snp'], map[,'feature'], combn, 2)) %>% 
            t %>% 
            as.data.frame
        ft <- graph_from_data_frame(ft, directed = FALSE)
        ft <- set_vertex_attr(ft, feature_type, index = match(map[,'snp'], V(ft)$name), 
                              map[,'feature'])
        
        n_features <- aggregate(feature ~ ., data=map, length)
        ft <- set_vertex_attr(ft, "n_features", value = 0)
        ft <- set_vertex_attr(ft, "n_features", 
                              index = match(n_features[,'snp'], V(ft)$name), 
                              n_features$feature)
    } else {
        stop("insufficient information to add feature information")
        return(NULL)
    }
    
    ft <- set_edge_attr(ft, "weight", value = 1)
    
    return(ft)
    
}

#' Map SNPs to genomic features.
#' 
#' @description Maps SNPs from a GWAS experiment to genomic features
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param flank A number with the flanking regions around features to be 
#' considered part of the feature i.e. SNPs mapped to them will be considered 
#' mapped to the feature.
#' @return A data.frame with two columns: one for the SNP and another for the
#' feature it has been mapped to.
#' @keywords internal
snp2feature <- function(gwas, features, flank = 0) {
    
    check_installed("IRanges", "snp2feature")
    
    # get map in appropriate format
    map <- gwas[["map"]]
    colnames(map) <- c("chr", "snp", "cm", "gpos", "allele1", "allele2")
    map[,'chr'] <- gsub("[Cc]hr", "", map[,'chr'])
    
    snp2feat <- by(map, map[,'chr'], function(snps) {
        
        # get features in the range defined by the snps
        grange <- paste(unique(snps[,'chr']), 
                        min(snps[,'gpos']) - flank, max(snps[,'gpos']) + flank, 
                        sep = ":")
        features <- subset(features, chromosome_name == unique(snps[,'chr']))
        
        if(nrow(features) == 0) {
            return()
        }
        
        # add a buffer before and after the feature
        features$start_position <- features$start_position - flank
        features$start_position[features$start_position < 0] <- 0
        features$end_position <- features$end_position + flank
        
        # convert to genomic ranges and check overlaps
        isnps <- with(snps, IRanges::IRanges(gpos, width=1, names=snp))
        ifeatures <- with(features, 
                       IRanges::IRanges(start_position, end_position, 
                                        names=feature))
        olaps <- IRanges::findOverlaps(isnps, ifeatures)
        s2f <- cbind(snps[S4Vectors::queryHits(olaps),], 
                     features[S4Vectors::subjectHits(olaps),])
        subset(s2f, select = c("snp", "feature"))
        
    })
    snp2feat <- do.call("rbind", snp2feat)
    
    return(snp2feat)
    
}