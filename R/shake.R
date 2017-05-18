shake <- function(gwas, net, ...) {

  # change gwas according to stuff
  # add some match correspondence between gwas and net
  # present the results better
  
  X <- as(gwas$genotypes, "numeric")
  Y <- gwas$fam$affected
  
  # remove redundant edges and self-edges
  net <- simplify(net)
  W <- get.adjacency(net, type="both", sparse = TRUE)
  settings <- getGinSettings(...)
  
  results <- runGin(X, Y, W, settings)
  
  return(results)
  
}