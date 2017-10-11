#' Find connected explanatory SNPs.
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being connected in an underlying network (Azencott et al., 2013).
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param net An igraph network that connects the SNPs.
#' @param ... Extra arguments for shake.
#' @return A copy of the SnpMatrix$map object, with the following additions:
#' \itemize{
#' \item{C: contains the univariate association score for every single SNP.}
#' \item{selected: logical vector indicating if the SNP was selected by shake or not.}
#' \item{module: integer with the number of the module the SNP belongs to.}
#' }
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & Borgwardt, K. M. (2013). Efficient network-guided multi-locus 
#' association mapping with graph cuts. Bioinformatics, 29(13), 171-179. \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom  igraph simplify as_adj
#' @importFrom methods as
#' @export
find_cones <- function(gwas, net, ...) {

  X <- as(gwas$genotypes, "numeric")
  Y <- gwas$fam$affected
  
  # remove redundant edges and self-edges
  net <- simplify(net)
  W <- as_adj(net, type="both", sparse = TRUE)
  
  # order according to order in map
  W <- W[gwas$map$snp.names, gwas$map$snp.names]
  
  settings <- get_shake_settings(...)
  
  gin <- run_shake(X, Y, W, settings)
  cat("eta =", gin$eta, "\nlambda =", gin$lambda, "\n")
  
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  map$C <- gin$scores
  map$selected <- as.logical(gin$indicator)
  
  map <- get_snp_modules(map, net)
  
  return(map)
  
}

#' Get shake settings.
#' 
#' @description Creates a list composed by all \code{shake} settings, with the values provided by the user, or the default ones if none is 
#' provided.
#' 
#' @param ... Any \code{shake} option.
#' @return A list of \code{shake} settings.
get_shake_settings <- function(...){
  
  settings <- list(...)
  
  # unsigned int
  if (! "nParameters" %in% names(settings))
    settings[["nParameters"]] = 10;
  
  # unsigned int
  if (! "folds" %in% names(settings))
    settings[["folds"]] = 10
  
  # bool, VectorXd, VectorXd
  if (! ("autoParameters" %in% names(settings) & "lambdas" %in% names(settings) & "etas" %in% names(settings)) ) {
    settings[["autoParameters"]] = TRUE
    settings[["lambdas"]] = numeric()
    settings[["lambdas"]] = rep(0, settings[["nParameters"]])
    settings[["etas"]] = numeric()
    settings[["etas"]] = rep(0, settings[["nParameters"]])
  }
  
  # unsigned int
  if (! "test_statistic" %in% names(settings))
    settings[["test_statistic"]] = 0
  
  # unsigned int
  if (! "gridsearch_depth" %in% names(settings))
    settings[["gridsearch_depth"]] = 1
  
  # unsigned int
  if (! "selection_criterion" %in% names(settings))
    settings[["selection_criterion"]] = 1
  
  # double
  if (! "seed" %in% names(settings))
    settings[["seed"]] = 0
  
  # double
  if (! "selection_ratio" %in% names(settings))
    settings[["selection_ratio"]] = 0.8
  
  # bool
  if (! "dump_intermediate_results" %in% names(settings))
    settings[["dump_intermediate_results"]] = TRUE
  
  # bool
  if (! "evaluateObjective" %in% names(settings))
    settings[["evaluateObjective"]] = FALSE
  
  return(settings);
}

