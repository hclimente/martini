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

#' Change encoding of dataset
#' 
#' @description Converts the encoding from additive to another one.
#' 
#' @param X Genotype matrix with SNPs encoded as 0 for major homozygous, 1 for
#' heterozygous and 2 for minor homozygous.
#' @param encoding Genetic model assumed: additive, recessive, dominant or
#' codominant.
#' @return A genotype matrix 
#' @examples 
#' X <- as(minigwas[["genotypes"]], "numeric")
#' martini:::encode_gwas(X, "recessive")
#' @keywords internal
encode_gwas <- function(X, encoding) {
  
  if (! encoding %in% c("additive", "recessive", "dominant", "codominant")) {
    stop("Invalid encoding.", call. = FALSE)
  }
  
  code <- data.frame(e = c('additive', 'recessive', 'dominant', 'codominant'),
                     AA = c(0,0,0,0),
                     AB = c(1,0,1,1),
                     BB = c(2,1,1,0))
  
  X[X == 0] <- code[code$e == encoding, 'AA']
  X[X == 1] <- code[code$e == encoding, 'AB']
  X[X == 2] <- code[code$e == encoding, 'BB']
  
  return(X)
  
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

#' Prepare covariates for \code{search_cones}
#' 
#' @description Prepares de covariates data.frame for the functions used in
#' \code{search_cones}, like \code{single_snp_association} or \code{score_folds}
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