#' Include LD information in the network.
#' 
#' @description Include LD information in the SNP network. Only work for human SNPs.
#' 
#' @param net A SNP network.
#' @param ld Data frame with linkage disequilibrium, like \code{get_ld} output.
#' @param method How to incorporate LD values into the network.
#' @return An SNP network where the edges weight 1 - LD, measured as Pearson correlation.
#' @importFrom igraph E %>% set_edge_attr delete_edges get.edgelist
#' @export
ldweight_edges <- function(net, ld, method = "inverse") {
  
  edges <- get.edgelist(net) %>% apply(1, sort) %>% t %>% apply(1, paste, collapse = "-")
  ld <- subset(ld, key %in% edges)
  
  snps <- strsplit(ld$key, "-") %>% unlist
  
  if (method == "inverse") {
    net <- set_edge_attr(net, "weight", index = E(net, P=snps), value = 1 / (1 + ld$r2))
  } else if (method == "subtraction") {
    net <- set_edge_attr(net, "weight", index = E(net, P=snps), value = 1 - ld$r2)
  }
  
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

#' Calculate LD in the datasets
#' 
#' @description Calculate LD as pairwise correlations in the GWAS dataset. Creates a sliding window for each SNP and calculates all the 
#' pairwise correlations inside that window.
#' 
#' @param gwas A GWAS experiment, in snpMatrix form.
#' @param window Window size.
#' @return An dataframe with a column key, with SNP ids separated by a "-" character, and an r2 column, 
#' containing the Pearson correlation.
#' @importFrom igraph %>%
#' @export
get_ld <- function(gwas, window = 5e4) {
  
  map <- gwas$map
  colnames(map) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  
  X <- as(gwas$genotypes, "numeric")
  
  # select only controls
  binary <- (unique(gwas$fam$affected) %>% length) == 2
  if (binary) {
    X <- X[gwas$fam$affected == 1, ]
  }
  
  ld <- by(map, map$chr, function(chr) {
    c <- unique(chr$chr)
    start <- min(chr$pos)
    end <- max(chr$pos)
    
    lapply(seq(start, end - window/2, window/2), function(x){
      
      mask <- (map$chr == c) & (map$pos >= x) & (map$pos <= x + window)
      
      if (sum(mask) < 2) {
        return(NULL)
      }
      
      r <- cor(X[,mask])
      r <- r[lower.tri(r)]
      
      snps <- as.character(map$snp[mask])
      
      keys <- lapply(1:(length(snps) - 1), function(x) {
        this <- snps[x]
        lapply(snps[(x+1):length(snps)], function(that) {
          sort(c(this, that)) %>% paste(collapse = "-")
        }) %>% unlist
      }) %>% unlist
      
      data.frame(key = keys, r2 = r^2)
      
    }) %>% invisible %>% do.call(rbind, .)
  }) %>% do.call(rbind, .) %>%
  unique
  
  ld$key <- as.character(ld$key)
  # remove cases where LD cannot be calculated, because one SNP presents no variance
  ld <- subset(ld, !is.na(r2))
  
  return(ld)
}