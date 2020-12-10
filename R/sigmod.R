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

  opts <- parse_scones_settings(c = 1, score, criterion, etas, lambdas)
  c <- single_snp_association(gwas, covars, opts[['score']])
  opts <- parse_scones_settings(c, score, criterion, etas, lambdas, TRUE)
  
  return(mincut.cv(gwas, net, covars, opts))
    
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
  
  return(mincut(gwas, net, covars, eta, lambda, score, sigmod = TRUE))
  
}
