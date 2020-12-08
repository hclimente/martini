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
sigmod.cv <- function(gwas, net, covars = data.frame(), score = "chi2", 
                      criterion = "consistency", etas = numeric(), 
                      lambdas = numeric()) {

  # get adjacency
  A <- get_A(gwas, net)
  
  # flip sign of lambdas
  opts <- parse_scones_settings(c = 1, score, criterion, etas, lambdas)
  c <- single_snp_association(gwas, covars, opts[['score']])
  opts <- parse_scones_settings(c = c, score, criterion, etas, lambdas)
  opts[['lambdas']] <- -opts[['lambdas']]
  
  return(mincut.cv(gwas, net, A, covars, opts))
    
}

#' @inherit scones
#' @references Liu, Y., Brossard, M., Roqueiro, D., Margaritte-Jeannin, P., 
#' Sarnowski, C., Bouzigon, E., Demenais, F. (2017). SigMod: an exact and
#' efficient method to identify a strongly interconnected disease-associated
#' module in a gene network. Bioinformatics, 33(10), 1536–1544. 
#' \url{https://doi.org/10.1093/bioinformatics/btx004}
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' sigmod(minigwas, gi, 10, 1)
#' @export
sigmod <- function(gwas, net, eta, lambda, covars = data.frame(), score = 'chi2') {
  
  A <- get_A(gwas, net)
  return(mincut(gwas, net, A, covars, eta, lambda, score))
  
}

#' Compute adjacency matrix
#' 
#' @template params_gwas
#' @template params_net
#' @return An adjacency matrix.
#' @importFrom igraph simplify as_adj
#' @importFrom Matrix diag rowSums
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

#' Parse \code{sigmod.cv} settings
#' 
#' @description Creates a list composed by all \code{sigmod.cv} settings, with 
#' the values provided by the user, or the default ones if none is provided.
#' @template params_gwas
#' @template params_covars
#' @template params_c
#' @template params_score
#' @template params_criterion
#' @template params_scones
#' @return A list of \code{sigmod.cv} settings.
#' @examples 
#' martini:::parse_sigmod_settings(minigwas, etas = c(1,2,3), lambdas = c(4,5,6))
#' martini:::parse_sigmod_settings(minigwas, c = c(1,10,100), score = "glm")
#' @keywords internal
parse_sigmod_settings <- function(gwas, covars = data.frame(),
                                  c = as.numeric(), ...) {
  
  if (!length(c)) {
    opts <- parse_scones_settings(c = 1, ...)
    c <- single_snp_association(gwas, covars, opts[['score']])
  }
  opts <- parse_scones_settings(c = c, ...)
  # flip sign of lambdas
  opts[['lambdas']] <- -opts[['lambdas']]
  
  return(opts)
}