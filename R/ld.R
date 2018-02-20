#' Include linkage disequilibrium information in the network.
#' 
#' @description Include linkage disequilibrium information in the SNP network. 
#' The weight of the edges will be lower the higher the linkage is.
#' 
#' @param net A SNP network.
#' @param ld A \code{dsCMatrix} or \code{dgCMatrix} containing linkage 
#' disequilibrium measures, like the output of \code{\link[snpStats]{ld}}.
#' @param method How to incorporate linkage-disequilibrium values into the 
#' network.
#' @return An copy of net where the edges weighted according to linkage 
#' disequilibrium.
#' @importFrom igraph E %>% set_edge_attr delete_edges get.edgelist
#' @importFrom Matrix summary
#' @importFrom stats cor
#' @examples 
#' ld <- snpStats::ld(minigwas$genotypes, depth = 2, stats = "R.squared")
#' # don't weight edges for which LD cannot be calculated
#' ld[is.na(ld)] <- 0
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' ldGi <- ldweight_edges(gi, ld)
#' @export
ldweight_edges <- function(net, ld, method = "inverse") {
    
    edges <- as.data.frame(get.edgelist(net))
    edges <- paste(edges[,1], edges[,2], sep = "-")
    ldDf <- Matrix::summary(ld)
    ldDf[['key']] <- paste(rownames(ld)[ldDf[['i']]], 
                           rownames(ld)[ldDf[['j']]], sep = "-")
    ldDf <- ldDf[ldDf[['key']] %in% edges, ]
    idx <- match(ldDf[['key']], edges)
    
    net <- switch(method,
        inverse = set_edge_attr(net, "weight", index = idx, 
                                value = 1 / (1 + ldDf[['x']])),
        subtraction = set_edge_attr(net, "weight", index = idx,
                                    value = 1 - ldDf[['x']]),
        sigmoid = set_edge_attr(net, "weight", index = idx, 
                                value = 1 / (1 + exp(10*(ldDf[['x']] - 0.5)))))

    if (any(is.na(E(net)$weight))) {
        stop("NA values as edge weights.")
    } else if (any(E(net)$weight < 0)) {
        stop("Edge weights cannot be negative.")
    }
    
    # remove edges with 0 weight
    zeroE <- E(net)[E(net)$weight == 0]
    net <- delete_edges(net, zeroE)
    
    return(net)
    
}