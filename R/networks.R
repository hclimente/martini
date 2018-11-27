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
#' @importFrom igraph %>% graph_from_data_frame simplify set_vertex_attr V
#' set_edge_attr
#' @importFrom utils combn
#' @examples 
#' get_GM_network(minigwas, snpMapping = minisnpMapping)
#' @export
get_GM_network <- function(gwas, organism = 9606, 
                           snpMapping = snp2gene(gwas, organism),
                           col_genes = c('snp','gene')) {
  
  snpMapping <- subset(snpMapping, select = col_genes)
  colnames(snpMapping) <- c('snp','gene')
  
  map <- gwas[["map"]]
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  map <- merge(map, snpMapping)
  map <- subset(map, select = c("snp","gene"))
  # in some cases the same position is linked to two different variants
  map <- unique(map)
  map[,'snp'] <- as.character(map[,'snp'])
  map[,'gene'] <- as.character(map[,'gene'])
  
  nSnps <- aggregate(snp ~ ., data=map, length)
  map <- subset(map, gene %in% nSnps[nSnps$snp > 1, 'gene'])
  
  gs <- get_GS_network(gwas)
  
  if (nrow(map) > 0) {
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

#' Map SNPs to genes.
#' 
#' @description Maps SNPs from a GWAS experiment to genes.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Organism: hsapiens represents human, athaliana for
#' Arabidopsis thaliana, etc.
#' @param flank A number with the flanking regions around genes to be considered
#' part of the gene i.e. SNPs mapped to them will be considered mapped to the
#' gene.
#' @return A data.frame with two columns: one for the SNP and another for the
#' gene it has been mapped to.
#' @keywords internal
snp2gene <- function(gwas, organism = 9606, flank = 0) {
  
  genes <- get_ensembl_genes(gwas, organism)
  snp2feature(gwas, genes, flank)
  
}

#' Get genomic position of the genes.
#' 
#' @description Get the genomic position of all the genes of an organism.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @return A data.frame with four columns, corresponding to the gene, the 
#' chromosome, the start and the end.
#' @keywords internal
get_ensembl_genes <- function(gwas, organism = 9606) {
    
  check_installed("biomaRt", "get_ensembl_genes")
  check_installed("httr", "get_ensembl_genes")
  check_installed("IRanges", "get_ensembl_genes")
  
  # get map in appropriate format
  map <- gwas[["map"]]
  colnames(map) <- c("chr", "snp", "cm", "gpos", "allele1", "allele2")
  map[,'chr'] <- gsub("[Cc]hr", "", map[,'chr'])
  
  # convert taxid to ensembl species name e.g. human databases are hsapiens_*
  urlTaxonomy <- "https://rest.ensembl.org/taxonomy"
  query <- paste0(urlTaxonomy,"/id/",organism,"?content-type=application/json")
  organism <- httr::GET(query)
  httr::stop_for_status(organism)
  organism <- httr::content(organism,type ="application/json",encoding="UTF-8")
  organism <- unlist(strsplit(organism$name, " "))
  organism <- tolower(paste0(substr(organism[1], 1,1), organism[2]))
  
  # create mart from ENSEMBL
  # consider vertebrates and plants
  ensembl <- tryCatch({
    datasetName <- paste0(organism, "_gene_ensembl")
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = datasetName)
  }, error = function(e){
    tryCatch({
      datasetName <- paste0(organism, "_eg_gene")
      biomaRt::useMart("plants_mart", host="plants.ensembl.org", 
               dataset=datasetName)  
    }, error = function(e) {
      stop(paste0("unable to find an appropriate database for ", organism, "."))
    })
  })
  
  genes <- biomaRt::getBM(
    attributes = c('ensembl_gene_id', 'external_gene_name', 
             'chromosome_name', 'start_position', 'end_position'),
    filters = 'biotype',
    values = list(biotype='protein_coding'), 
    mart = ensembl)
  noName <- genes[,'external_gene_name'] == "" | 
    is.na(genes[,'external_gene_name'])
  genes$feature <- ifelse(noName, genes[,'ensembl_gene_id'], genes[,'external_gene_name'])
  genes <- subset(genes, select = c('feature', 'chromosome_name', 
                                    'start_position', 'end_position'))
  
  return(genes)
  
}

#' Get protein-protein interactions.
#' 
#' @description Get all protein-protein interactions for an organism from
#' BioGRID.
#' 
#' @param organism Organism: human represents human, arabidopsis for Arabidopsis
#' thaliana, etc.
#' @return A data.frame with two columns with pairs of interacting proteins.
#' @examples 
#' # download dog interactions
#' martini:::get_ppi(9615)
#' @keywords internal
get_ppi <- function(organism = 9606) {
  
  check_installed("httr", "get_ppi")
  
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
    colnames(p) <- c("geneA","geneB")
    return(p)
  })
  
  ppi <- do.call("rbind", ppi)
  ppi <- unique(ppi)
  
  return(ppi)
  
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
        
        # get genes in the range defined by the snps
        grange <- paste(unique(snps[,'chr']), 
                        min(snps[,'gpos']) - flank, max(snps[,'gpos']) + flank, 
                        sep = ":")
        features <- subset(features, chromosome_name == unique(snps[,'chr']))
        
        if(nrow(features) == 0) {
            return()
        }
        
        # add a buffer before and after the gene
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