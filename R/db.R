#' Map SNPs to Ensembl genes.
#' 
#' @description Maps SNPs from a GWAS experiment to genes.
#' @template params_gwas
#' @template params_organism
#' @param flank A number with the flanking regions around genes to be considered
#' part of the gene i.e. SNPs mapped to them will be considered mapped to the
#' gene.
#' @return A data.frame with two columns: one for the SNP and another for the
#' gene it has been mapped to.
#' @keywords internal
snp2ensembl <- function(gwas, organism = 9606, flank = 0) {
  
  check_installed(c("IRanges", "httr", "biomaRt"), "snp2ensembl")
  
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
#' @template params_organism
#' @return A data.frame with two columns with pairs of interacting proteins.
#' @examples 
#' # download dog interactions
#' \dontrun{martini:::get_gxg_biogrid(9615)}
#' @keywords internal
get_gxg_biogrid <- function(organism = 9606) {
  
  check_installed(c("httr"), "get_gxg_biogrid")
  
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
#' 
#' @description Get all protein-protein interactions for an organism from
#' STRING. It uses a score cut-off of 400.
#' @template params_organism
#' @return A data.frame with two columns with pairs of interacting proteins.
#' @importFrom igraph as_edgelist
#' @examples 
#' # download frog interactions
#' \dontrun{martini:::get_gxg_string(8364)}
#' @keywords internal
get_gxg_string <- function(organism = 9606) {
  
  check_installed(c("STRINGdb", "httr", "biomaRt"), "get_gxg_string")
  
  string_db <- STRINGdb::STRINGdb$new(version = '11.5', species = organism,
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
#' @inherit get_gxg_biogrid
#' @keywords internal
mget_gxg_biogrid <- memoise::memoise(get_gxg_biogrid)

#' Memoised version of get_gxg_stringdb
#' @inherit get_gxg_string
#' @keywords internal
mget_gxg_string <- memoise::memoise(get_gxg_string)

#' Get gene interactions
#' 
#' @description Wrapper for the different functions to get gene-gene 
#' interactions. Supports cached results.
#' @param db String containing the database to obtain the gene-gene interactions
#' from. Possible values: 'biogrid', 'string'.
#' @template params_organism
#' @template params_flush
#' @return A data.frame with two columns with pairs of interacting proteins.
#' @keywords internal
get_gxg <- function(db, organism, flush) {
  
  check_installed(c("biomaRt", "httr", "memoise", "STRINGdb"), "get_gxg")
  
  if (db == 'biogrid') {
    f <- mget_gxg_biogrid
  } else if (db == 'string') {
    f <- mget_gxg_string
  } else {
    stop(paste('unknown gene interaction database', db))
  }
  
  if (memoise::has_cache(f)(organism)) {
    if (flush) {
      message('cache flushed!')
      memoise::drop_cache(f)(organism)
    } else {
      warning('using cache. Use flush = TRUE to get new gene interactions.')
    }
  }
  
  return(f(organism))
  
}
