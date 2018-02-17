#' Calculate an empirical p-value for \code{find_cones} results.
#' 
#' @description Do a permutation-based test to assess the statistical
#' significance of each of the modules obtained through shake. For a module of 
#' size k, k interconnected SNPs are picked \code{nperm} times, and their joint 
#' association score is calculated to come up with an estimation of the 
#' distribution.
#' 
#' @param cones Results from \code{find_cones}.
#' @param net The same SNP network provided to \code{find_cones}.
#' @param nperm Integer with the name of permutations.
#' @return An empirical p-value for each of the SNP modules. Please, note that 
#' the minimum possible p-value from an empirical distribution is set to 
#' \code{1/(nperm+1)}. Modules composed by a single SNP will not have an 
#' empirical p-value.
#' @importFrom stats ecdf
#' @importFrom igraph vcount random_walk
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' cones <- search_cones(minigwas, gi)
#' test_cones_modules(cones, gi, 100)
#' @export
test_cones_modules <- function(cones, net, nperm = 100000) {
    
    numSNPs <- vcount(net)
    selected <- subset(cones, selected)
    
    # calculate one ecdf for each module size
    moduleSizes <- unique(table(selected$module))
    names(moduleSizes) <- moduleSizes
    
    ecdfs <- lapply(moduleSizes, function(n) {
        if (n > 1){
            sampled_C <- lapply(seq_len(nperm), function(i){
                v <- sample(seq_len(numSNPs), 1)
                snpmodule_i <- random_walk(net, v, n)
                cones_module <- cones[cones$snp %in% snpmodule_i, ]
                sum(cones_module$c)
            })
            sampled_C <- do.call("c", sampled_C)
            
            ecdf(sampled_C)
        }
    })
    
    modules <- by(selected, selected$module, function(k){
        n <- nrow(k)
        i <- unique(k$module)
        
        if (n > 1){
            snpmodule_C <- sum(k$c)
            p <- 1 - ecdfs[[as.character(n)]](snpmodule_C)
        } else {
            p <- NA
        }
        
        data.frame(k = i, ncomponents = n, p = p)
    })
    
    modules <- do.call("rbind", modules)
    # minimum p-vaue = 1 / (nperm + 1)
    modules$p[modules$p == 0] = 1 / (nperm + 1)
    
    return(modules)
    
}