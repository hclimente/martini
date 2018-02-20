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
#' @importFrom igraph V induced_subgraph vertex_attr %>%
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' martini:::subnet(gi, "gene", "A")
#' martini:::subnet(gi, "name", c("1A1", "1A3"))
#' @keywords internal
subnet <- function(net, attr, values, affirmative = TRUE) {
  vertices <- V(net)[(vertex_attr(net, attr) %in% values) == affirmative]
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
  V(net)[(vertex_attr(net, attr) %in% values) == affirmative]
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
  
  X[X == 0] <- code$AA[code$e == encoding]
  X[X == 1] <- code$AB[code$e == encoding]
  X[X == 2] <- code$BB[code$e == encoding]
  
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