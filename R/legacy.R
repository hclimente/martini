#' Find connected explanatory SNPs.
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being
#'  connected in an underlying network (Azencott et al., 2013).
#' @template params_gwas
#' @template params_net
#' @param encoding SNP encoding (unused argument).
#' @param sigmod Boolean. If TRUE, use the Sigmod variant of SConES, meant to 
#' prioritize tightly connected clusters of SNPs.
#' @template params_covars
#' @param associationScore Association score to measure association between 
#' genotype and phenotype.
#' @param modelScore String with the function to measure the quality of a split.
#' @template params_etas
#' @template params_lambdas
#' @template return_cones
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & 
#' Borgwardt, K. M. (2013). Efficient network-guided multi-locus 
#' association mapping with graph cuts. Bioinformatics, 29(13), 171-179. 
#' \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph simplify as_adj
#' @importFrom Matrix diag rowSums
#' @importFrom methods as
#' @examples
#' \dontrun{gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' search_cones(minigwas, gi)
#' search_cones(minigwas, gi, encoding = "recessive")
#' search_cones(minigwas, gi, associationScore = "skat")}
#' @export
search_cones <- function(gwas, net, encoding = "additive", sigmod = FALSE,
                         covars = data.frame(), 
                         associationScore = c("chi2", "glm"), 
                         modelScore = c("stability", "bic", "aic", "aicc", 
                                        "global_clustering", "local_clustering"), 
                         etas = numeric(), lambdas = numeric()) {
  
  .Deprecated('scones.cv')
  
  if (!missing("encoding"))
    warning("Argument encoding ignored, will be set to additive.")
  
  associationScore <- match.arg(associationScore)
  modelScore <- match.arg(modelScore)
  if (modelScore == 'consistency') modelScore <- 'stability'
  
  args <- list(gwas=gwas, net=net, covars=covars, score=associationScore, 
               criterion=modelScore, etas=etas, lambdas=lambdas)
  if(sigmod) cones <- do.call(sigmod.cv, args)
  else cones <- do.call(scones.cv, args)
  
  return(cones)
  
}
