#' Subgraph of vertices with an attribute
#' 
#' @description Returns a subgraph matching some condition.
#' 
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}
#' @return A subgraph containing only the vertices with attribute equal to any
#' of the values in \code{values}.
#' @importFrom igraph V induced_subgraph vertex_attr %>%
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' martini:::subnet(gi, "gene", "A")
#' martini:::subnet(gi, "name", c("1A1", "1A3"))
subnet <- function(net, attr, values) {
  vertices <- V(net)[vertex_attr(net, attr) %in% values]
  induced_subgraph(net, vertices)
}

#' Vertices with an attribute
#' 
#' @description Returns the nodes matching some condition.
#' 
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}
#' @return The vertices with attribute equal to any of the values in
#' \code{values}.
#' @importFrom igraph V vertex_attr
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' martini:::subvert(gi, "gene", "A")
#' martini:::subvert(gi, "name", c("1A1", "1A3"))
subvert <- function(net, attr, values) {
  V(net)[vertex_attr(net, attr) %in% values]
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
#' X <- as(minigwas$genotypes, "numeric")
#' martini:::encode_gwas(X, "recessive")
encode_gwas <- function(X, encoding) {
  
  if (encoding == "additive") {
    AA <- 0; AB <- 1; BB <- 2
  } else if (encoding == "recessive") {
    AA <- 0; AB <- 0; BB <- 1
  } else if (encoding == "dominant") {
    AA <- 0; AB <- 1; BB <- 1
  } else if (encoding == "codominant") {
    AA <- 0; AB <- 1; BB <- 0
  } else {
    stop("Invalid encoding.", call. = FALSE)
  }
  
  X[X == 0] <- AA
  X[X == 1] <- AB
  X[X == 2] <- BB
  
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
is_coherent <- function(gwas) {
  
  mapSelfCoherence <- by(gwas$map, gwas$map[,1], function(chr) {
    ascendingChr <- chr[order(chr[,4]),]
    descendingChr <- chr[order(chr[,4], decreasing = TRUE),]
    return(any(chr != ascendingChr) & any(chr != descendingChr))
  })
  
  if (any(mapSelfCoherence)) {
    stop("map is not ordered by genomic position.")
  }
  
  if (any(gwas$map[,2] != colnames(gwas$genotypes))) {
    stop("map and genotype SNP order differ.")
  }
  
  return(TRUE)
  
}