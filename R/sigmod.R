#' @inherit scones.cv
#' @template params_scones
#' @references Liu, Y., Brossard, M., Roqueiro, D., Margaritte-Jeannin, P., 
#' Sarnowski, C., Bouzigon, E., Demenais, F. (2017). SigMod: an exact and
#' efficient method to identify a strongly interconnected disease-associated
#' module in a gene network. Bioinformatics, 33(10), 1536–1544. 
#' \url{https://doi.org/10.1093/bioinformatics/btx004}
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' sigmod.cv(minigwas, gi)
#' sigmod.cv(minigwas, gi, score = "glm")
#' @export
sigmod.cv <- function(gwas, net, covars = data.frame(), ...) {

  map <- sanitize_map(gwas)
  
  # get adjacency
  net <- simplify(net)
  A <- as_adj(net, type="both", sparse = TRUE, attr = "weight")
  A <- W[map[['snp']], map[['snp']]]
  diag(A) <- -rowSums(A)
  
  # flip sign of lambdas
  opts <- parse_scones_settings(c = 1, ...)
  c <- single_snp_association(gwas, covars, opts[['score']])
  opts <- parse_scones_settings(c = c, ...)
  opts[['lambdas']] <- -opts[['lambdas']]
  
  return(mincut.cv(gwas, net, A, covars, opts))
    
}

#' @inherit scones
#' @export
sigmod <- function(gwas, net, eta, lambda, score = 'chi2', covars = data.frame()) {
  
  return(scones(gwas, net, eta, -lambda, score, covars))
  
}