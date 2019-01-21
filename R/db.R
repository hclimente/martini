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