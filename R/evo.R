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
#' \item{selected: logical vector indicating if the SNP was selected by evo or not.}
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
#' @param cones Results from \code{evo}.
#' @param net The same SNP network provided to \code{evo}.
#' @return A list with the clusters of selected SNPs.
#' @importFrom igraph induced_subgraph components
cluster_snps <- function(cones, net) {
  
  selected <- subset(cones, selected)
  subnet <- induced_subgraph(net, selected$snp)
  
  clusters <- components(subnet)
  cones$cluster <- NA
  cones$cluster[cones$selected] <- clusters$membership[order(match(names(clusters$membership), cones$snp))]
  
  return(cones)
}

#' Get evo settings.
#' 
#' @description Creates a list composed by all \code{evo} settings, with the values provided by the user, or the default ones if none is 
#' provided.
#' @param associationScore Association score to measure association between genotype and phenotype. Possible values: chi2 (default), skat, 
#' trend.
#' @param modelScore Model selection criterion Possible values: consistency, bic (default), aic, aicc, mbic.
#' @param encoding SNP encoding. Possible values: additive (default), resessive, dominant, codominant.
#' @param etas Numeric vector with the etas to explore in the grid search. If ommited, it's automatically created based on the association
#' scores.
#' @param lambdas Numeric vector with the lambdas to explore in the grid search. If ommited, it's automatically created based on the 
#' association scores.
#' @param debug Display additional information. Possible values: TRUE, FALSE (default).
#' @return A list of \code{evo} settings.
get_evo_settings <- function(associationScore = "chi2", modelScore = "bic", encoding = "additive", etas = numeric(), lambdas = numeric(), 
                             debug = FALSE){
  
  settings <- list()
  
  # unsigned int
  if (associationScore == "skat") {
    settings$associationScore = 0
  } else if (associationScore == "chi2") {
    settings$associationScore = 1
  } else if (associationScore == "trend") {
    settings$associationScore = 2
  } else {
    stop(paste("Error: invalid associationScore", associationScore))
  }
  
  # unsigned int
  if (modelScore == "consistency") {
    settings$modelScore = 0
  } else if (modelScore == "bic") {
    settings$modelScore = 1
  } else if (modelScore == "aic") {
    settings$modelScore = 2
  } else if (modelScore == "aicc") {
    settings$modelScore = 3
  } else if (modelScore == "mbic") {
    settings$modelScore = 4
  } else {
    stop(paste("Error: invalid modelScore", modelScore))
  }
  
  # unsigned int
  if (encoding == "additive") {
    settings$encoding = 0
  } else if (encoding == "recessive") {
    settings$encoding = 1
  } else if (encoding == "dominant") {
    settings$encoding = 2
  } else if (encoding == "codominant") {
    settings$encoding = 3
  } else {
    stop(paste("Error: invalid encoding", encoding))
  }
  
  # bool
  if (! is.logical(debug)) {
    stop("Error: debug must be logical.")
  } else {
    settings$debug = debug;
  }
  
  # VectorXd
  if (!is.numeric(etas)){
    stop("Error: etas must be numeric")
  } else {
    settings$etas = etas
  }
  
  # VectorXd
  if (!is.numeric(lambdas)){
    stop("Error: lambdas must be numeric")
  } else {
    settings$lambdas = lambdas
  }
  
  return(settings);
}
