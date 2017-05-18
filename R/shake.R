shake <- function(gwas, net, ...) {

  # add some match correspondence between gwas and net
  
  X <- as(gwas$genotypes, "numeric")
  Y <- gwas$fam$affected
  
  # remove redundant edges and self-edges
  net <- simplify(net)
  W <- get.adjacency(net, type="both", sparse = TRUE)
  settings <- getGinSettings(...)
  
  gin <- runGin(X, Y, W, settings)
  
  gwas$map$ginscore <- gin$scores
  gwas$map$ginpicked <- as.logical(gin$indicator)
  gwas$gin <- list(lambda = gin$lambda, eta = gin$eta)
  
  return(gwas)
  
}