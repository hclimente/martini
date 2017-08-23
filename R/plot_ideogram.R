#' Ideogram of shake results.
#' 
#' @description Create a circular ideogram of the \code{shake} results using the circlize package (Gu et al., 2014).
#' 
#' @param map Output from \code{shake}.
#' @param net The SNP network provided to \code{shake}.
#' @return A circular ideogram, including the manhattan plot, and the interactions between the selected SNPs.
#' @references Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014). circlize Implements and enhances circular visualization in R. Bioinformatics (Oxford, England), 30(19), 2811-2. \url{https://doi.org/10.1093/bioinformatics/btu393}
#' @importFrom circlize circos.initializeWithIdeogram circos.genomicTrackPlotRegion circos.genomicPoints circos.genomicLink circos.clear
#' @importFrom igraph get.data.frame delete_vertices
#' @export
plot_ideogram <- function(map, net){
  
  circos.initializeWithIdeogram()
  
  bed <- map2bed(map)
  bed$C <- map$C
  bed$selected <- map$selected
  # order to put the selected snps in fron in the plot
  bed <- bed[with(bed, order(selected)),]
  
  circos.genomicTrackPlotRegion(bed, 
                                ylim = c(0, 1.1 * max(bed$C)), 
                                panel.fun = function(region, value, ...) {
                                  # color according to selection/non-selection
                                  col = ifelse(value[[2]], "orange", "gray70")
                                  circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
                                }, track.height = 0.3)
  
  # create links
  # remove unselected SNPs from the graph
  df <- delete_vertices(net, as.character(map$snp[! map$selected])) %>%
    get.data.frame %>%
    subset(select = c("from", "to"))
  
  region1 <- merge(df, map, by.x = "from", by.y = "snp", sort = F) %>%
    subset(select = c("chr", "from", "cm", "pos")) %>%
    map2bed
  region2 <- merge(df, map, by.x = "to", by.y = "snp", sort = F) %>%
    subset(select = c("chr", "to", "cm", "pos")) %>%
    map2bed
  
  circos.genomicLink(region1, region2)
  
  circos.clear()
  
}

#' Converts a MAP data.frame to a BED data.frame
#' 
#' @description Takes a map file and:
#'  \itemize{
#' \item{column 1: Used as the chromosome column in the BED file..}
#' \item{column 4: Used as start and end in the BED data.frame (as we work with SNPs).}
#' }
#' 
#' @param map A MAP data.frame.
#' @return A BED data.frame.
map2bed <- function(map) {
  
  bed <- map[,c(1,4)]
  colnames(bed) <- c("chr", "start")
  bed$chr <- paste0("chr", bed$chr)
  bed$chr <- ifelse(bed$chr == "chr23", "chrX", bed$chr)
  bed$end <- bed$start
  
  return(bed)
}