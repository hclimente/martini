#' Include LD information in the network.
#' 
#' @description Include LD information in the SNP network. Only work for human SNPs.
#' 
#' @param net A SNP network.
#' @return An SNP network where the edges weight 1 - LD, measured as Pearson correlation.
#' @importFrom igraph E %>% set_edge_attr delete_edges get.edgelist
#' @export
ldweight_edges <- function(net, ld = get_ld(net)) {
  
  edges <- get.edgelist(net) %>% apply(1, sort) %>% t %>% apply(1, paste, collapse = "_")
  ld <- subset(ld, key %in% edges)
  
  snps <- strsplit(ld$key, "_") %>% unlist
  
  net <- set_edge_attr(net, "weight", index = E(net, P=snps), value = 1 - ld$r2)
  
  if (any(is.na(E(net)$weight))) {
    stop("NA values as edge weights.")
  } else if (any(E(net)$weight < 0)) {
    stop("Pearson coefficients cannot be > 1.")
  }
  
  # remove edges with 0 weight
  zeroE <- E(net)[E(net)$weight == 0]
  net <- delete_edges(net, zeroE)
  
  return(net)
  
}

get_ld <- function(net = NULL, gwas = NULL, organism = 9606, db = "1000GENOMES:phase_3:CEU") {
  
  if (!is.null(gwas)) {
    ld <- get_ld_from_gwas(gwas, organism, db)
  } else if (!is.null(net)) {
    ld <- get_ld_from_net(net, organism, db)
  } else {
    stop("You must specify either gwas or net.")
  }
  
  return(ld)
}

#' Retrieve LD information from Ensembl
#' 
#' @description Retrieve LD information from Ensembl in the format required by ldweight_edges.
#' 
#' @param gwas A GWAS experiment, in snpMatrix form.
#' @return An dataframe with a column key, with SNP ids separated by an underscore, and an r2 column, 
#' containing the Pearson correlation.
#' @importFrom httr GET content stop_for_status content_type
#' @importFrom igraph %>%
get_ld_from_gwas <- function(gwas, organism, db) {
  
  baseUrl <- "https://rest.ensembl.org/ld"
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  
  ld <- by(map, map$chr, function(C) {
    
    start <- min(C$pos)
    end <- max(C$pos)
    chr <- unique(C$chr)
    
    lapply(seq(start, end, 1000000), function(chunk) {
      
      c <- subset(C, pos >= chunk & pos <= chunk + 1000000)
      
      if (nrow(c) > 1) {
        
        region <- paste0(chr, ":", min(c$pos), "..", max(c$pos) )
        q <- paste(baseUrl, organism, "region", region, db, sep = "/")
        r <- GET(q, content_type("application/json"))
        stop_for_status(r)
        r <- content(r, type="application/json", encoding="UTF-8", simplifyDataFrame = T)
        
        r <- subset(r, variation1 %in% c$snp & variation2 %in% c$snp)
        
        if (nrow(r) > 0) {
          r$key <- apply(r[,c("variation1","variation2")], 1, sort) %>% t %>% apply(1, paste, collapse = "_")
          r <- subset(r, select = c("key", "r2"))
          return(r)
        }
      }
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)

  ld <- unique(ld)
  ld$r2 <- as.numeric(ld$r2)
  
  return(ld)
  
}

#' Retrieve LD information from Ensembl
#' 
#' @description Retrieve LD information from Ensembl in the format required by ldweight_edges.
#' 
#' @param net A SNP network.
#' @return An dataframe with a column key, with SNP ids separated by an underscore, and an r2 column, 
#' containing the Pearson correlation.
#' @importFrom igraph V %>% get.edgelist
#' @importFrom httr GET content stop_for_status content_type
get_ld_from_net <- function(net, organism, db) {
  
  baseUrl <- "https://rest.ensembl.org/ld"
  edges <- get.edgelist(net) %>% apply(1, sort) %>% t %>% apply(1, paste, collapse = "_")
  
  ld <- lapply(names(V(net)), function(snp) {
    q <- paste(baseUrl, organism, snp, db, sep = "/")
    r <- GET(q, content_type("application/json"))
    stop_for_status(r)
    r <- content(r, type="application/json", encoding="UTF-8", simplifyDataFrame = T)
    
    if (length(r) > 0) {
      r$key <- apply(r[,c("variation1","variation2")], 1, sort) %>% t %>% apply(1, paste, collapse = "_")
      r <- subset(r, key %in% edges, select = c("key", "r2"))
      return(r)
    }
  }) %>% do.call(rbind, .)
  
  ld <- unique(ld)
  ld$r2 <- as.numeric(ld$r2)
  
  return(ld)
  
}