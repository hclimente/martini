#' Make pseudo SnpMatrix object
#' 
#' @description Wrap design matrix and outcome vector into a pseudo SnpMatrix
#' object.
#' @template params_X
#' @template params_y
#' @template params_featnames
#' @template params_net
wrap_Xy <- function(X, y, featnames, net) {
  
  featnames_2 <- intersect(featnames, names(V(net)))
  
  gwas <- list('genotypes' = X[, featnames %in% featnames_2], 
               'fam' = data.frame(fam = 0,
                                  affected = y),
               'map' = data.frame(chr = 0, 
                                  snp = featnames[featnames %in% featnames_2], 
                                  cm = 0, 
                                  pos = 0, 
                                  allele.1 = 0, 
                                  allele.2 = 0))
  
  net <- induced_subgraph(net, featnames_2)
  net <- set_edge_attr(net, "weight", value=1)
  
  return(list(gwas = gwas, net = net))
  
}

#' Find connected explanatory features
#' 
#' @description Finds the features maximally associated with a phenotype while 
#' being connected in an underlying network. Select the hyperparameters by
#' cross-validation.
#' @template params_X
#' @template params_y
#' @template params_featnames
#' @template params_net
#' @template return_cones
#' @examples 
#' X <- as(minigwas[['genotypes']], 'numeric')
#' X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' scones.cv_(X, minigwas[['fam']]$affected, minigwas[['map']]$snp, gi)
#' @export
scones.cv_ <- function(X, y, featnames, net) {
  
    i <- wrap_Xy(X, y, featnames, net)
    scones.cv(i[['gwas']], i[['net']], score = 'r2')
  
}

#' Find connected explanatory features
#' 
#' @description Finds the features maximally associated with a phenotype while 
#' being connected in an underlying network.
#' @template params_X
#' @template params_y
#' @template params_featnames
#' @template params_net
#' @template params_eta
#' @template params_lambda
#' @template return_cones
#' @examples 
#' X <- as(minigwas[['genotypes']], 'numeric')
#' X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' scones_(X, minigwas[['fam']]$affected, minigwas[['map']]$snp, gi, 10, 1)
#' @export
scones_ <- function(X, y, featnames, net, eta, lambda) {
  
    i <- wrap_Xy(X, y, featnames, net)
    scones(i[['gwas']], i[['net']], eta, lambda, score = 'r2')
  
}

#' @inherit scones.cv_
#' @examples 
#' X <- as(minigwas[['genotypes']], 'numeric')
#' X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' sigmod.cv_(X, minigwas[['fam']]$affected, minigwas[['map']]$snp, gi)
#' @export
sigmod.cv_ <- function(X, y, featnames, net) {
  
  i <- wrap_Xy(X, y, featnames, net)
  sigmod.cv(i[['gwas']], i[['net']], score = 'r2')
  
}

#' @inherit scones_
#' @examples 
#' X <- as(minigwas[['genotypes']], 'numeric')
#' X <- X + matrix(rnorm(2500, sd = 0.1), nrow(X), ncol(X))
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' sigmod_(X, minigwas[['fam']]$affected, minigwas[['map']]$snp, gi, 10, 1)
#' @export
sigmod_ <- function(X, y, featnames, net, eta, lambda) {
  
  i <- wrap_Xy(X, y, featnames, net)
  sigmod(i[['gwas']], i[['net']], eta, lambda, score = 'r2')
  
}
