#' Subgraph of vertices with an attribute
#' 
#' @description Returns a subgraph matching some condition.
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}.
#' @param affirmative Logical. States if a condition must be its affirmation 
#' (e.g. all nodes with gene name "X"), or its negation (all nodes not with gene
#' name "X").
#' @return A subgraph containing only the vertices with attribute equal to any
#' of the values in \code{values}.
#' @importFrom igraph induced_subgraph
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' martini:::subnet(gi, "gene", "A")
#' martini:::subnet(gi, "name", c("1A1", "1A3"))
#' @keywords internal
subnet <- function(net, attr, values, affirmative = TRUE) {
  vertices <- subvert(net, attr, values, affirmative)
  induced_subgraph(net, vertices)
}

#' Vertices with an attribute
#' 
#' @description Returns the nodes matching some condition.
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}
#' @param affirmative Logical. States if a condition must be its affirmation 
#' (e.g. all nodes with gene name "X"), or its negation (all nodes not with gene
#' name "X").
#' @return The vertices with attribute equal to any of the values in
#' \code{values}.
#' @importFrom igraph V vertex_attr
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' martini:::subvert(gi, "gene", "A")
#' martini:::subvert(gi, "name", c("1A1", "1A3"))
subvert <- function(net, attr, values, affirmative = TRUE) {
  
  n_attr <- vertex_attr(net, attr)
  select <- (n_attr %in% values) == affirmative
  
  if (affirmative) {
    select <- !is.na(n_attr) & select
  }
  
  V(net)[select]
}

#' Check package is installed
#' 
#' @description Checks if a package is installed, launches an error if it is
#' not.
#' 
#' @param pkgs Character vector with the names of the packages.
#' @param fn Function calling the check.
#' @return The package is loaded into the namespace.
#' @examples 
#' martini:::check_installed(c("martini"))
#' \dontrun{martini:::check_installed("martinid")}
#' @keywords internal
check_installed <- function(pkgs, fn = "This function") {
  installed <- unlist(lapply(pkgs, requireNamespace, quietly = TRUE))
  if (!all(installed)) {
    stop(paste0(fn, " requires the following packages to be installed:\n", 
                paste(pkgs[!installed], collapse = '\n')), call. = FALSE)
  }
}

#' Check inner coherence of GWAS dataset
#' 
#' @description Checks that the different data structures have the SNPs in the 
#' same order.
#' @template params_gwas
#' @return TRUE if the GWAS dataset is coherent. Else, raises an error.
#' @examples 
#' martini:::is_coherent(minigwas)
#' @keywords internal
is_coherent <- function(gwas) {
  
  mapSelfCoherence <- by(gwas[["map"]], gwas[["map"]][,1], function(chr) {
    ascendingChr <- chr[order(chr[,4]),]
    descendingChr <- chr[order(chr[,4], decreasing = TRUE),]
    return(any(chr != ascendingChr) & any(chr != descendingChr))
  })
  
  if (any(mapSelfCoherence)) {
    stop("map is not ordered by genomic position.")
  }
  
  if (any(gwas[["map"]][,2] != colnames(gwas[["genotypes"]]))) {
    stop("map and genotype SNP order differ.")
  }
  
  return(TRUE)
  
}

#' Prepare covariates for \code{scones}
#' 
#' @description Prepares de covariates data.frame for the functions used in
#' \code{scones}, like \code{single_snp_association} or \code{score_folds}
#' .
#' @template params_gwas
#' @template params_covars
#' @return The covars data.frame, with the rows in the same order as gwas. 
#' @keywords internal
arrange_covars <- function(gwas, covars) {
  
  if (any(dim(covars))) {
    covars <- covars[match(row.names(gwas[['genotypes']]), covars[['sample']]),]
    covars <- subset(covars, select = -sample)
  }
  
  return(covars)

}

#' Permute samples
#' 
#' @description Compute a permutation of the samples of a snpMatrix object. 
#' Useful to make sure that the folds are not stratified by phenotype.
#' @template params_gwas
#' @keywords internal
permute_snpMatrix <- function(gwas) {
  
  n <- nrow(gwas[['genotypes']])
  perm <- sample(n)
  
  gwas[['genotypes']] <- gwas[['genotypes']][perm,]
  gwas[['fam']] <- gwas[['fam']][perm,]
  
  return(gwas)
  
}

#' Subsample snpMatrix
#' 
#' @description Compute a permutation of the samples of a snpMatrix object. 
#' Useful to make sure that the folds are not stratified by phenotype.
#' @template params_gwas
#' @param samples Vector (logical or numeric) containing the samples to select.
#' @keywords internal
subset_snpMatrix <- function(gwas, samples) {
  
  gwas[['genotypes']] <- gwas[['genotypes']][samples,]
  gwas[['fam']] <- gwas[['fam']][samples,]
  
  return(gwas)
  
}


#' Check snpMapping
#' 
#' @description Check that snpMapping is a proper data.frame.
#' 
#' @param snpMapping data.frame containing the correspondence between SNPs and
#' genes.
#' @param col_genes Length 2 character vector containing the colnames containing
#' the SNP and the gene ids, respectively.
#' @keywords internal
sanitize_snpMapping <- function(snpMapping, col_genes) {
  
  if(!is.null(snpMapping) && nrow(snpMapping)) {
    snpMapping <- subset(snpMapping, select = col_genes)
  } else {
    warning("no mappings between SNPs and genes were provided.")
    snpMapping <- data.frame(matrix(ncol = 2, nrow = 0))
  }
  colnames(snpMapping) <- c('snp','gene')
  
  return(snpMapping)
  
}

#' Check map
#' 
#' @description Check that map is a proper data.frame.
#' @template params_gwas
#' @keywords internal
sanitize_map <- function(gwas) {
  
  map <- gwas[["map"]]
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  
  return(map)
  
}

#' Tax id to ensembl species name
#' @description Converts taxid to ensembl species name e.g. human databases are 
#' hsapiens_*
#' @template params_organism
#' @keywords internal
organism_id2name <- function(id) {
  
  check_installed(c("httr"), "organism_id2name")
  
  urlTaxonomy <- "http://rest.ensembl.org/taxonomy"
  query <- paste0(urlTaxonomy,"/id/",id,"?content-type=application/json")
  name <- httr::GET(query)
  httr::stop_for_status(name)
  name <- httr::content(name,type ="application/json",encoding="UTF-8")
  name <- unlist(strsplit(name$name, " "))
  name <- tolower(paste0(substr(name[1], 1,1), name[2]))
  
  return(name)
  
}

#' Open a biomaRt connection
#' @description Opens a biomaRt connection for the relevant species.
#' @param organism String containing the ensembl species name (e.g. hsapiens 
#' for human)
#' @keywords internal
connect_biomart <- function(organism) {
  
  check_installed(c("biomaRt"), "connect_biomart")
  
  # create mart from ENSEMBL
  # consider vertebrates and plants
  conn <- tryCatch({
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
  
  return(conn)
  
}

#' Compute Laplacian matrix
#' 
#' @template params_gwas
#' @template params_net
#' @return A Laplacian matrix.
#' @importFrom igraph simplify as_adj
#' @importFrom Matrix diag rowSums
#' @keywords internal
get_adjacency <- function(gwas, net) {
  
  map <- sanitize_map(gwas)
  
  # remove redundant edges and self-edges in network and sort
  net <- simplify(net)
  L <- as_adj(net, type="both", sparse = TRUE, attr = "weight")
  L <- L[map[['snp']], map[['snp']]]
  
  return(L)
  
}

#' Converts a MAP data.frame to a BED data.frame
#'
#' @description Takes a map file and:
#' \itemize{
#' \item{column 1: Used as the chromosome column in the BED file.}
#' \item{column 4: Used as start and end in the BED data.frame (as we work with
#' SNPs).}
#' }
#' @template params_gwas
#' @return A BED data.frame.
gwas2bed <- function(gwas) {
  
  map <- sanitize_map(gwas)
  bed <- subset(map, select = c("chr", "pos"))
  colnames(bed) <- c("chr", "start")
  bed$chr <- paste0("chr", bed$chr)
  bed$chr <- ifelse(bed$chr == "chr23", "chrX", bed$chr)
  bed$end <- bed$start
  
  return(bed)
  
}
