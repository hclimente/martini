#' Subgraph of vertices with an attribute
#' 
#' @description Returns a subgraph matching some condition.
#' 
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
#' 
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
#' @param pkg Name of the package.
#' @param fn Function calling the check.
#' @return The package is loaded into the namespace.
#' @examples 
#' martini:::check_installed("martini")
#' \dontrun{martini:::check_installed("martinid")}
#' @keywords internal
check_installed <- function(pkg, fn = "this function") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(pkg, "needed for", fn, "to work. Please install it."),
         call. = FALSE)
  }
}

#' Check inner coherence of GWAS dataset
#' 
#' @description Checks that the different data structures have the SNPs in the 
#' same order.
#' @param gwas A GWAS experiment.
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
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param net An igraph network that connects the SNPs.
#' @keywords internal
arrange_covars <- function(gwas, covars) {
  
  if (ncol(covars)) {
    covars <- covars[match(row.names(gwas[['genotypes']]), covars[['sample']]), ]
    covars <- subset(covars, select = -sample)
  }
  
  return(covars)

}

#' Permute samples
#' 
#' @description Compute a permutation of the samples of a snpMatrix object. 
#' Useful to make sure that the folds are not stratified by phenotype.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @keywords internal
permute_snpMatrix <- function(gwas) {
  
  n <- nrow(gwas[['genotypes']])
  perm <- sample(n)
  
  gwas[['genotypes']] <- gwas[['genotypes']][perm,]
  gwas[['fam']] <- gwas[['fam']][perm,]
  
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
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @keywords internal
sanitize_map <- function(gwas) {
  
  map <- gwas[["map"]]
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  
  return(map)
  
}