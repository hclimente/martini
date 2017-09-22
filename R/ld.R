#' Include LD information in the network.
#' 
#' @description Include LD information in the SNP network. Only work for human SNPs.
#' 
#' @param net A SNP network.
#' @return An SNP network where the edges weight 1 - LD, measured as Pearson correlation.
#' @importFrom igraph E %>% set_edge_attr delete_edges
#' @export
ldweight_edges <- function(net, ld = get_ld(net)) {
  
  snps <- strsplit(ld$key, "_") %>% unlist
  
  net <- set_edge_attr(net, "weight", index = E(net, P = snps), value = 1 - ld$r2)
  
  if (any(E(net)$weight < 0)) {
    stop("Pearson coefficients cannot be > 1.")
  }
  
  # remove edges with 0 weight
  zeroE <- E(net)[E(net)$weight == 0]
  net <- delete_edges(net, zeroE)
  
  return(net)
  
}

#' Retrieve LD information from Ensembl
#' 
#' @description Retrieve LD information from Ensembl in the format required by ldweight_edges.
#' 
#' @param net A SNP network.
#' @return An dataframe with a column key, with SNP ids separated by an underscore, and an r2 column, 
#' containing the Pearson correlation.
#' @importFrom igraph V %>% get.edgelist
#' @importFrom httr GET content stop_for_status
#' @export
get_ld <- function(net) {
  
  baseUrl <- paste0("https://rest.ensembl.org/ld/human/")
  outputType <- "/1000GENOMES:phase_3:KHV?content-type=application/json"
  edges <- get.edgelist(net) %>% apply(1, sort) %>% t %>% apply(1, paste, collapse = "_")
  
  ld <- lapply(names(V(net)), function(snp) {
    r <- GET(paste0(baseUrl, snp, outputType))
    stop_for_status(r)
    r <- content(r, type="application/json", encoding="UTF-8", simplifyDataFrame = T)
    
    if (length(r) > 0) {
      r$key <- apply(r[,c("variation1","variation2")], 1, sort) %>% t %>% apply(1, paste, collapse = "_")
      r <- subset(r, key %in% edges, select = c("key", "r2"))
      return(r)
    }
  })
  ld <- do.call(rbind, ld)
  ld <- unique(ld)
  ld$r2 <- as.numeric(ld$r2)
  
  return(ld)
  
}