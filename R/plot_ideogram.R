#' Ideogram of shake results.
#' 
#' @description Create a circular ideogram of the \code{shake} results using the circlize package (Gu et al., 2014).
#' 
#' @param results Results from \code{shake}.
#' @param net The same SNP network provided to \code{shake}.
#' @return A circular ideogram, including the manhattan plot, and the interactions between the selected SNPs.
#' @references Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014). circlize Implements and enhances circular visualization in R. Bioinformatics (Oxford, England), 30(19), 2811–2. \url{https://doi.org/10.1093/bioinformatics/btu393}
#' @export
plot_ideogram <- function(results, net){
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("circlize needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  circos.initializeWithIdeogram()
  
  bed <- map2bed(results)
  bed$C <- results$C
  bed$selected <- results$selected
  # order to put the selected snps in fron in the plot
  bed <- bed[with(bed, order(selected)),]
  
  circos.genomicTrackPlotRegion(bed, ylim = c(0, 1.1 * max(bed$C)), panel.fun = function(region, value, ...) {
    # color according to selection/non-selection
    col = ifelse(value[[2]], "orange", "gray70")
    circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
  }, track.height = 0.3)
  
  # create links
  # remove unselected SNPs from the graph
  net <- delete_vertices(net, results$snp[! results$selected])
  df <- get.data.frame(net)
  
  region1 <- merge(df, results, by.x = "from", by.y = "snp", sort = F)
  region1 <- map2bed(region1)
  region2 <- merge(df, results, by.x = "to", by.y = "snp", sort = F)
  region2 <- map2bed(region2)
  
  circos.genomicLink(region1, region2)
  
  circos.clear()
  
}

map2bed <- function(map) {
  
  bed <- subset(map, select = c("chr", "pos"))
  colnames(bed) <- c("chr", "start")
  bed$chr <- paste0("chr", bed$chr)
  bed$chr <- ifelse(bed$chr == "chr23", "chrX", bed$chr)
  bed$end <- bed$start
  
  return(bed)
}