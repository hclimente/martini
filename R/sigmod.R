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

  # get adjacency
  A <- get_A(gwas, net)
  
  # flip sign of lambdas
  opts <- parse_scones_settings(c = 1, ...)
  c <- single_snp_association(gwas, covars, opts[['score']])
  opts <- parse_scones_settings(c = c, ...)
  opts[['lambdas']] <- -opts[['lambdas']]
  
  return(mincut.cv(gwas, net, A, covars, opts))
    
}

#' @inherit scones
#' @export
sigmod <- function(gwas, net, eta, lambda, covars = data.frame(), score = 'chi2') {
  
  A <- get_A(gwas, net)
  return(mincut(gwas, net, A, covars, eta, lambda, score))
  
}

#' @keywords internal
get_A <- function(gwas, net) {
  
  map <- sanitize_map(gwas)
  
  # remove redundant edges and self-edges in network and sort
  net <- simplify(net)
  A <- as_adj(net, type="both", sparse = TRUE, attr = "weight")
  A <- A[map[['snp']], map[['snp']]]
  diag(A) <- -rowSums(A)
  
  return(A)
  
}

#' @keywords internal
parse_sigmod_settings <- function(gwas, covars, ...) {
  
  opts <- parse_scones_settings(c = 1, ...)
  c <- single_snp_association(gwas, covars, opts[['score']])
  opts <- parse_scones_settings(c = c, ...)
  # flip sign of lambdas
  opts[['lambdas']] <- -opts[['lambdas']]
  
  return(opts)
}