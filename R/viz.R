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
#' @importFrom igraph as_ids
#' @export
plot_ideogram <- function(gwas, net, covars = data.frame(), genome = "hg19") {
  
  check_installed(c("circlize", "IRanges"), "plot_ideogram")
  
  circlize::circos.initializeWithIdeogram(species = genome)
  
  bed <- gwas2bed(gwas)
  bed[['snp']] <- sanitize_map(gwas)[['snp']]
  bed[['c']] <- snp_test(gwas, covars, 'chi2')
  bed[['selected']] <- bed[['snp']] %in% names(V(net))
  # bring the selected snps to the front
  bed <- bed[with(bed, order(selected)),]
  
  circlize::circos.genomicTrackPlotRegion(
    bed,
    ylim = c(0, 1.1 * max(bed$c, na.rm = TRUE)),
    panel.fun = function(region, value, ...) {
      # color according to selection/non-selection
      col = ifelse(value[[3]], "#e84545", "#bad7df")
      circlize::circos.genomicPoints(region, value, col = col,
                                     cex = 0.5, pch = 16)
    }, 
    track.height = 0.3)
  
  # create edges
  edges <- as_ids(E(net))
  edges <- do.call(rbind, strsplit(edges, '|', fixed = TRUE))
  edges <- as.data.frame(edges)
  
  ## get positional information and remove self interactions
  edges <- merge(edges, bed, by.x = 'V1', by.y = 'snp')
  edges <- merge(edges, bed, by.x = 'V2', by.y = 'snp')
  edges <- edges[edges[['chr.x']] != edges[['chr.y']],]
  
  # group into regions
  edges <- cbind(group_snps(edges, "chr.x", "start.x", 50000),
                 group_snps(edges, "chr.y", "start.y", 50000))
  edges <- unique(edges)
  
  x1 <- edges[,c(1,2,3)]
  x2 <- edges[,c(4,5,6)]
  
  circlize::circos.genomicLink(x1, x2, col=sample(1:5, nrow(x1), replace=TRUE))
  circlize::circos.clear()
  
}

#' Groups nearby SNPs
#'
#' @description Groups SNPs closer than a specifiec threshold of distance.
#' @param bed data.frame containing at least two properties (chromosome and 
#' position) of a set of SNPs.
#' @param chr_col Name of the column containing the SNP chromosome.
#' @param pos_col Name of the column containing the SNP position.
#' @param threshold Maximum distance to group two SNPs group.
#' @return A data.frame in bed format, with the same dimensions as the original,
#' but with the groups.
#' @keywords internal
group_snps <- function(bed, chr_col, pos_col, threshold) {
  
  check_installed("IRanges", "group_snps")
  
  # make row numbers to reorder at the end
  bed[['id']]  <- 1:nrow(bed)
  
  ir_bed <- by(bed, bed[,chr_col], function(chr) {
    ir_chr <- IRanges::IRanges(start = chr[,pos_col], width = threshold)
    reduced_chr <- IRanges::reduce(ir_chr)
    ol <- as.matrix(IRanges::findOverlaps(ir_chr, reduced_chr))
    ir_chr[ol[,'queryHits'],] <- reduced_chr[ol[,'subjectHits'],]
    ir_chr <- as.data.frame(ir_chr)[,c('start','end')]
    ir_chr <- data.frame(chr_range   = chr[,chr_col],
                         start_range = ir_chr[['start']],
                         end_range   = ir_chr[['end']])
    unique(cbind(ir_chr, chr[,c(chr_col, pos_col)]))
          
  }) %>% do.call(rbind, .)
  
  bed <- merge(bed, ir_bed, by = c(chr_col, pos_col))
  bed[order(bed[['id']]), c('chr_range','start_range','end_range')]
  
}