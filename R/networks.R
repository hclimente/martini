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
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Tax ID of the studied organism. Required if snpMapping is not
#' provided. The default is 9606 (human).
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
#' @importFrom igraph %>% graph_from_data_frame simplify set_vertex_attr V
#' set_edge_attr
#' @importFrom utils combn
#' @examples 
#' get_GM_network(minigwas, snpMapping = minisnpMapping)
#' @export
get_GM_network <- function(gwas, organism = 9606, 
                           snpMapping = snp2gene(gwas, organism),
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
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Tax ID of the studied organism. Required if snpMapping is not
#' provided. The default is 9606 (human).
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
#' @param flush Remove cached results? Boolean value.
#' @return An igraph network of the GI network of the SNPs.
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., &
#' Borgwardt, K. M. (2013). Efficient network-guided multi-locus association
#' mapping with graph cuts. Bioinformatics, 29(13), 171-179.
#' \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph graph_from_data_frame simplify set_edge_attr
#' @examples 
#' get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' @export
get_GI_network <- function(gwas, organism = 9606,
                           snpMapping = snp2gene(gwas, organism), 
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

#' Map SNPs to genes.
#' 
#' @description Maps SNPs from a GWAS experiment to genes.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Organism: 9606 represents human, etc.
#' @param flank A number with the flanking regions around genes to be considered
#' part of the gene i.e. SNPs mapped to them will be considered mapped to the
#' gene.
#' @return A data.frame with two columns: one for the SNP and another for the
#' gene it has been mapped to.
#' @keywords internal
snp2gene <- function(gwas, organism = 9606, flank = 0) {
  
  check_installed("IRanges", "snp2gene")
  
  # get map in appropriate format
  map <- sanitize_map(gwas)
  map[,'chr'] <- gsub("[Cc]hr", "", map[,'chr'])
  
  organism <- organism_id2name(organism)
  ensembl <- connect_biomart(organism)
  
  snp2gene <- by(map, map[,'chr'], function(snps) {
    
    # get genes in the range defined by the snps
    grange <- paste(unique(snps[,'chr']), 
                    min(snps[,'pos']) - flank, max(snps[,'pos']) + flank, 
                    sep = ":")
    genes <- biomaRt::getBM(
                   attributes = c("ensembl_gene_id","external_gene_name",
                                  "start_position","end_position"),
                   filters = c("chromosomal_region", "biotype"),
                   values = list(chromosomal_region=grange,
                                 biotype="protein_coding"), 
                   mart = ensembl)
    
    if(nrow(genes) == 0) {
      return()
    }
    
    # add a buffer before and after the gene
    genes$start_position <- genes$start_position - flank
    genes$start_position[genes$start_position < 0] <- 0
    genes$end_position <- genes$end_position + flank
    
    # convert to genomic ranges and check overlaps
    isnps <- with(snps, IRanges::IRanges(pos, width=1, names=snp))
    igenes <- with(genes, 
                   IRanges::IRanges(start_position, end_position, 
                                    names=ensembl_gene_id))
    olaps <- IRanges::findOverlaps(isnps, igenes)
    s2g <- cbind(snps[S4Vectors::queryHits(olaps),], 
                 genes[S4Vectors::subjectHits(olaps),])
    hasName <- s2g[,'external_gene_name'] == "" | 
               is.na(s2g[,'external_gene_name'])
    s2g$gene <- ifelse(hasName, s2g[,'ensembl_gene_id'], 
                                s2g[,'external_gene_name'])
    subset(s2g, select = c("snp", "gene"))
    
  })
  snp2gene <- do.call("rbind", snp2gene)
  
  return(snp2gene)
  
}

#' Get BioGRID protein-protein interactions.
#' 
#' @description Get all protein-protein interactions for an organism from
#' BioGRID.
#' 
#' @param organism Organism: human represents human, arabidopsis for Arabidopsis
#' thaliana, etc.
#' @return A data.frame with two columns with pairs of interacting proteins.
#' @examples 
#' # download dog interactions
#' martini:::get_gxg_biogrid(9615)
#' @keywords internal
get_gxg_biogrid <- function(organism = 9606) {
  
  check_installed("httr", "get_gxg_biogrid")
  
  # construct query: all interactions in the requested organism
  baseUrl <- "http://webservice.thebiogrid.org/interactions/?"
  query <- paste(baseUrl,
                 "accessKey=912520d90dba1547c01abf31a18bcc94",
                 "interSpeciesExcluded=true", 
                 "includeHeader=true", 
                 paste0("taxId=", organism), sep = "&")
  
  # number of results
  N <- httr::GET(paste(query, "format=count", sep = "&"))
  httr::stop_for_status(N)
  N <- as.numeric(
             httr::content(N, type="text/csv", encoding="UTF-8", col_types="i"))
  
  # retrieve results in batches
  ppi <- lapply(seq(1, N, 10000), function(i){
    q <- paste(query, paste0("start=", i), sep = "&")
    req <- httr::GET(q)
    httr::stop_for_status(req)
    
    # parse results
    biogrid <- httr::content(req, type="text/tab-separated-values", 
                             encoding="UTF-8", 
                             col_types="cccccccccccccccccccccccc")
    p <- subset(biogrid, select=c("Official Symbol Interactor A",
                                  "Official Symbol Interactor B"))
    colnames(p) <- c("gene1","gene2")
    return(p)
  })
  
  ppi <- do.call("rbind", ppi)
  ppi <- unique(ppi)
  
  return(ppi)
  
}

#' Get STRING protein-protein interactions.
#' @description Get all protein-protein interactions for an organism from
#' STRING. It uses a score cut-off of 400.
#' @param organism Organism: 9606 represents human, etc.
#' @return A data.frame with two columns with pairs of interacting proteins.
#' @importFrom igraph as_edgelist
#' @examples 
#' # download frog interactions
#' martini:::get_gxg_string(8364)
#' @keywords internal
get_gxg_string <- function(organism = 9606) {
  
  check_installed("STRINGdb", "get_gxg_string")
  
  string_db <- STRINGdb::STRINGdb$new(version = '11', species = organism,
                                      score_threshold = 400)
  
  # download ppi
  ppi <- string_db$get_graph()
  ppi <- as_edgelist(ppi)
  ppi <- as.data.frame(ppi)
  
  ppi[['V1']] <- gsub(paste0(organism, '.'), '', ppi[['V1']])
  ppi[['V2']] <- gsub(paste0(organism, '.'), '', ppi[['V2']])
  
  # convert protein ids to gene ids
  organism <- organism_id2name(organism)
  ensembl <- connect_biomart(organism)
  
  genes <- unique(c(ppi[['V1']], ppi[['V2']]))
  ensp2hgnc <- biomaRt::getBM(filters = "ensembl_peptide_id", 
        attributes = c("ensembl_peptide_id", "hgnc_symbol"),
        values = genes, mart = ensembl)
  ensp2hgnc <- ensp2hgnc[ensp2hgnc[["hgnc_symbol"]] != "",]
  
  ppi <- merge(ppi, ensp2hgnc, by.x = 'V1', by.y = "ensembl_peptide_id")
  ppi <- merge(ppi, ensp2hgnc, by.x = 'V2', by.y = "ensembl_peptide_id")
  ppi <- ppi[,c('hgnc_symbol.x', 'hgnc_symbol.y')]
  ppi <- unique(ppi)
  colnames(ppi) <- c("gene1","gene2")
  
  return(ppi)
}

#' Memoised version of get_gxg_biogrid
#' 
#' @importFrom memoise memoise
#' @keywords internal
mget_gxg_biogrid <- memoise(get_gxg_biogrid)

#' Memoised version of get_gxg_stringdb
#' 
#' @importFrom memoise memoise
#' @keywords internal
mget_gxg_string <- memoise(get_gxg_string)

#' Get gene interactions
#' 
#' @description Wrapper for the different functions to get gene-gene 
#' interactions. Supports cached results.
#' @param db String containing the database to obtain the gene-gene interactions
#' from. Possible values: 'biogrid', 'string'.
#' @param organism Organism: 9606 represents human, etc.
#' @param flush Remove cached results? Boolean value.
#' @importFrom memoise memoise drop_cache has_cache
#' @keywords internal
get_gxg <- function(db, organism, flush) {
  
  if (db == 'biogrid') {
    f <- mget_gxg_biogrid
  } else if (db == 'string') {
    f <- mget_gxg_string
  } else {
    stop(paste('unknown gene interaction database', db))
  }
  
  if (has_cache(f)(organism)) {
    if (flush) {
      message('cache flushed!')
      drop_cache(f)(organism)
    } else {
      warning('using cache. Use flush = TRUE to get new gene interactions.')
    }
  }
  
  return(f(organism))
  
}