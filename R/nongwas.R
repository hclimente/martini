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
#' @export
scones.cv_ <- function(X, y, featnames, net) {
  
    i <- wrap_Xy(X, y, featnames, net)
    scones.cv(i[['gwas']], i[['net']], score = 'ttest')
  
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
#' @export
scones_ <- function(X, y, featnames, net, eta, lambda) {
  
    i <- wrap_Xy(X, y, featnames, net)
    scones(i[['gwas']], i[['net']], eta, lambda, score = 'ttest')
  
}

#' @inherit scones.cv_
#' @export
sigmod.cv_ <- function(X, y, featnames, net) {
  
  i <- wrap_Xy(X, y, featnames, net)
  sigmod.cv(i[['gwas']], i[['net']], score = 'ttest')
  
}

#' @inherit scones_
#' @export
sigmod_ <- function(X, y, featnames, net, eta, lambda) {
  
  i <- wrap_Xy(X, y, featnames, net)
  sigmod(i[['gwas']], i[['net']], eta, lambda, score = 'ttest')
  
}