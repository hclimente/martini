#' Get list of protein-protein interactions from IntAct.
#' 
#' @description Get all protein-protein interactions for an organism.
#' 
#' @param organism Organism: human represents human, arabidopsis for Arabidopsis thaliana, etc.
#' @return A dataframe with two columns with pairs of interacting proteins.
#' @export
get_ppi <- function(organism = "arabidopsis") {
  
  if (!requireNamespace("httr", quietly = TRUE))
    stop("httr needed for this function to work. Please install it.", call. = FALSE)
  
  # construct query: all interactions in the requested organism
  baseUrl <- "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/query/"
  query <- paste0(baseUrl, paste0("identifier:taxidA:", organism, "%20AND%20taxidB:", organism))
  
  # number of results
  N <- GET(paste0(query, "?format=count"))
  N <- as.numeric(content(N, encoding = "UTF-8"))
  
  # retrieve results in batches
  ppi <- lapply(seq(1, N, 1000), function(i){
    q <- paste0(query, "?format=tab25&firstResult=", i, "&maxResult=", min(i+999, N))
    req <- GET(q)
    
    # parse results
    intact <- content(req, type = "text/tab-separated-values", encoding = "UTF-8", col_names = F)
    
    # get genenames from aliases columns (5 and 6)
    p <- apply(subset(intact, select = c("X5", "X6")), 2, function(alias) {
      alias <- strsplit(alias, split = '|', fixed = T)
      alias <- lapply(alias, function(x) grep("uniprotkb:.*\\(gene name\\)", x, value = T) )
      alias <- lapply(alias, function(x) gsub("uniprotkb:", "", x, fixed = T) )
      alias <- lapply(alias, function(x) gsub("(gene name)", "", x, fixed = T) )
      alias
    })
    
    # keep values only if both have a genename
    p <- mapply(function(x,y){
      if(length(x) > 0 & length(y) > 0) {
        data.frame(geneA = x, geneB = y)
      }
    }, p$X5, p$X6)
    do.call("rbind", p)
  })
  
  ppi <- do.call("rbind", ppi)
  ppi <- unique(ppi)
  
  return(ppi)

}