#' @inherit scones.cv
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
sigmod.cv <- function(gwas, net, covars = data.frame(),
                      score = c("chi2", "glm", "r2"), 
                      criterion = c("stability", "bic", "aic", "aicc", 
                                    "global_clustering", "local_clustering"), 
                      etas = numeric(), lambdas = numeric(),
                      family = c("binomial", "poisson", "gaussian", "gamma"), 
                      link = c("logit", "log", "identity", "inverse"),
                      max_prop_snp = 0.5) {

  score <- match.arg(score)
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  link <- match.arg(link)
  c <- snp_test(gwas, covars, score, family, link)
  grid <- get_grid(c = c, etas, lambdas)
  
  return(mincut.cv(gwas, net, covars, grid[['etas']], grid[['lambdas']], 
                   criterion, score, TRUE, family, link, max_prop_snp))
  
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
sigmod <- function(gwas, net, eta, lambda, covars = data.frame(), 
                   score = c("chi2", "glm", "r2"), 
                   family = c("binomial", "poisson", "gaussian", "gamma"), 
                   link = c("logit", "log", "identity", "inverse")) {
  
  score <- match.arg(score)
  family <- match.arg(family)
  link <- match.arg(link)
  
  return(mincut(gwas, net, covars, eta, lambda, score, TRUE, family, link))
  
}
