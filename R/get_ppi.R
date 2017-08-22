#' Get protein-protein interactions.
#' 
#' @description Get all protein-protein interactions for an organism from BioGRID
#' 
#' @param organism Organism: human represents human, arabidopsis for Arabidopsis thaliana, etc.
#' @return A dataframe with two columns with pairs of interacting proteins.
#' @export
get_ppi <- function(organism = 9606) {
  
  # construct query: all interactions in the requested organism
  baseUrl <- "http://webservice.thebiogrid.org/interactions/?"
  query <- paste(baseUrl,
                 "accessKey=912520d90dba1547c01abf31a18bcc94",
                 "interSpeciesExcluded=true", 
                 "includeHeader=true", 
                 paste0("taxId=", organism), sep = "&")
  
  # number of results
  N <- httr::GET(paste(query, "format=count", sep = "&"))
  N <- as.numeric(httr::content(N, type="text/csv", encoding="UTF-8", col_types="i"))
  
  # retrieve results in batches
  ppi <- lapply(seq(1, N, 10000), function(i){
    q <- paste(query, paste0("start=", i), sep = "&")
    req <- httr::GET(q)
    
    # parse results
    biogrid <- httr::content(req, type="text/tab-separated-values", encoding="UTF-8", col_types="cccccccccccccccccccccccc")
    p <- subset(biogrid, select=c("Official Symbol Interactor A", "Official Symbol Interactor B"))
    colnames(p) <- c("geneA","geneB")
    return(p)
  })
  
  ppi <- do.call("rbind", ppi)
  ppi <- unique(ppi)
  
  return(ppi)

}