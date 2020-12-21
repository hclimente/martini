#' Ideogram of SConES results.
#'
#' @description Create a circular ideogram of the a network results using the
#' circlize package (Gu et al., 2014).
#' @template params_gwas
#' @template params_net
#' @template params_covars
#' @template params_genome
#' @return A circular ideogram, including the manhattan plot, and the
#' interactions between the selected SNPs.
#' @references Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014).
#' circlize Implements and enhances circular visualization in R. Bioinformatics
#' (Oxford, England), 30(19), 2811-2. 
#' \url{https://doi.org/10.1093/bioinformatics/btu393}
#' @export
plot_ideogram <- function(gwas, net, covars = data.frame(), genome = "hg19") {
  
  check_installed(c("circlize"), "plot_ideogram")
  
  circlize::circos.initializeWithIdeogram(species = genome)
  
  modules <- get_snp_modules(gwas, net)
  modules[['selected']] <- modules[['snp']] %in% names(V(net))

  bed <- gwas2bed(gwas)
  bed[['c']] <- snp_test(gwas, covars, 'chi2')
  bed[['selected']] <- modules[['snp']] %in% names(V(net))
  # bring the selected snps to the front
  bed <- bed[with(bed, order(selected)),]
  
  circlize::circos.genomicTrackPlotRegion(
    bed,
    ylim = c(0, 1.1 * max(bed$c, na.rm = T)),
    panel.fun = function(region, value, ...) {
      # color according to selection/non-selection
      col = ifelse(value[[2]], "#e84545", "#bad7df")
      circlize::circos.genomicPoints(region, value, col = col,
                                     cex = 0.5, pch = 16)
    }, 
    track.height = 0.3)
  
  # create links
  selected <- subset(modules, selected)
  
  regions <- by(selected, selected[,c("chr","module")], function(k) {
    data.frame(chr = paste0("chr", unique(k$chr)),
               module = unique(k$module),
               start = min(k$pos),
               end = max(k$pos))
  }) %>% do.call(rbind, .)
  
  links <- merge(regions, regions, by = "module")
  
  region1 <- subset(links, chr.x != chr.y, select = c("chr.x","start.x","end.x"))
  region2 <- subset(links, chr.x != chr.y, select = c("chr.y","start.y","end.y"))
  
  circlize::circos.genomicLink(region1, region2, 
                               col = sample(1:5, nrow(region1), replace = TRUE))
  
  circlize::circos.clear()
  
}
