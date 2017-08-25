#' Find connected explanatory SNPs.
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being connected in an underlying network (Azencott et al., 2013).
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param net An igraph network that connects the SNPs.
#' @param ... Extra arguments for \code{\link{get_evo_settings}}.
#' @return A copy of the SnpMatrix$map object, with the following additions:
#' \itemize{
#' \item{C: contains the univariate association score for every single SNP.}
#' \item{selected: logical vector indicating if the SNP was selected by shake or not.}
#' \item{cluster: integer with the number of the cluster the SNP belongs to.}
#' }
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus 
#' association mapping with graph cuts. Bioinformatics, 29(13), 171-179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph simplify as_adj
#' @export
search_cones <- function(gwas, net, ...) {

  X <- as(gwas$genotypes, "numeric")
  Y <- gwas$fam$affected
  
  # remove redundant edges and self-edges
  net <- igraph::simplify(net)
  W <- igraph::as_adj(net, type="both", sparse = TRUE)
  
  # order according to order in map
  W <- W[gwas$map$snp.names, gwas$map$snp.names]
  
  settings <- get_evo_settings(...)
  
  test <- evo(X, Y, W, settings)
  cat("eta =", test$eta, "\nlambda =", test$lambda, "\n")
  
  cones <- gwas$map
  colnames(cones) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  cones$c <- test$c
  cones$selected <- as.logical(test$selected)
  
  cones <- cluster_snps(cones, net)
  
  return(cones)
  
}

#' Return groups of interconnected SNPs.
#' 
#' @description Find clusters composed by interconnected SNPs.
#' 
#' @param map Results from \code{shake}.
#' @param net The same SNP network provided to \code{shake}.
#' @return A list with the clusters of selected SNPs.
#' @importFrom igraph induced_subgraph components
cluster_snps <- function(map, net) {
  
  selected <- subset(map, selected)
  subnet <- induced_subgraph(net, selected$snp)
  
  clusters <- components(subnet)
  map$cluster <- NA
  map$cluster[map$selected] <- clusters$membership[order(match(names(clusters$membership), map$snp))]
  
  return(map)
}

#' Get shake settings.
#' 
#' @description Creates a list composed by all \code{shake} settings, with the values provided by the user, or the default ones if none is 
#' provided.
#' @param associationScore Association score to measure association between genotype and phenotype. Possible values: skat (default), chi2, 
#' trend.
#' @param modelScore Model selection criterion Possible values: consistency, bic (default), aic, aicc, mbic.
#' @param encoding SNP encoding. Possible values: additive (default), resessive, dominant, codominant.
#' @param debug Display additional information. Possible values: TRUE, FALSE (default).
#' @return A list of \code{shake} settings.
get_evo_settings <- function(...){
  
  settings <- list(...)
  
  # unsigned int
  if (! "associationScore" %in% names(settings)) {
    settings[["associationScore"]] = 0
  } else if (settings[["associationScore"]] == "skat") {
    settings[["associationScore"]] = 0
  } else if (settings[["associationScore"]] == "chi2") {
    settings[["associationScore"]] = 1
  } else if (settings[["associationScore"]] == "trend") {
    settings[["associationScore"]] = 2
  } else {
    stop(paste("Error: invalid associationScore", settings[["associationScore"]]))
  }
  
  # unsigned int
  if (! "modelScore" %in% names(settings)) {
    settings[["modelScore"]] = 1
  } else if (settings[["modelScore"]] == "consistency") {
    settings[["modelScore"]] = 0
  } else if (settings[["modelScore"]] == "bic") {
    settings[["modelScore"]] = 1
  } else if (settings[["modelScore"]] == "aic") {
    settings[["modelScore"]] = 2
  } else if (settings[["modelScore"]] == "aicc") {
    settings[["modelScore"]] = 3
  } else if (settings[["modelScore"]] == "mbic") {
    settings[["modelScore"]] = 4
  } else {
    stop(paste("Error: invalid modelScore", settings[["modelScore"]]))
  }
  
  # unsigned int
  if (! "encoding" %in% names(settings)) {
    settings[["encoding"]] = 0;
  } else if (settings[["encoding"]] == "additive") {
    settings[["encoding"]] = 0
  } else if (settings[["encoding"]] == "recessive") {
    settings[["encoding"]] = 1
  } else if (settings[["encoding"]] == "dominant") {
    settings[["encoding"]] = 2
  } else if (settings[["encoding"]] == "codominant") {
    settings[["encoding"]] = 3
  } else {
    stop(paste("Error: invalid encoding", settings[["encoding"]]))
  }
  
  # bool
  if (! "debug" %in% names(settings)) {
    settings[["debug"]] = FALSE;
  } else if (! is.logical(settings[["debug"]])) {
    stop("Error: debug must be logical.")
  }
  
  return(settings);
}
