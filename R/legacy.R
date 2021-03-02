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
#' @importFrom igraph V
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
  if(sigmod) subnet <- do.call(sigmod.cv, args)
  else subnet <- do.call(scones.cv, args)
  
  cones <- get_snp_modules(gwas, subnet)
  cones[['c']] <- snp_test(gwas, covars, associationScore)
  cones[['selected']] <- cones[['snp']] %in% names(V(subnet))
  
  return(cones)
  
}

#' Return groups of interconnected SNPs.
#' 
#' @description Find modules composed by interconnected SNPs.
#' 
#' @template params_gwas
#' @template params_net
#' @return A list with the modules of selected SNPs.
#' @importFrom igraph induced_subgraph components
#' @examples
#' \dontrun{
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' cones <- search_cones(minigwas, gi)
#' martini:::get_snp_modules(cones, gi)
#' }
#' @keywords internal
get_snp_modules <- function(gwas, net) {
  
  # compute components
  modules <- components(net)
  modules <- as.data.frame(modules['membership'])
  colnames(modules) <- "module"
  modules['snp'] <- rownames(modules)
  
  # annotate on the map file
  map <- sanitize_map(gwas) 
  modules <- merge(map, modules, all.x = TRUE)
  map <- modules[match(map[,'snp'], modules[,'snp']),]
  
  return(map)
  
}
