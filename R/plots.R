#' Ideogram of a SNP module.
#' 
#' @description Create a ideogram of a SNP module from \code{SConES} results using the Gviz package (Hahne and Ivanek, 2016).
#' 
#' @param cones Results from \code{SConES}.
#' @param k Id of the module to plot.
#' @param genome Abbreviations of the genome to use: hg19 for human (default),  mm10 for mouse, etc. Argument to be passed to 
#' \code{\link{IdeogramTrack}}.
#' @return An ideogram per chromosome showing the selected SNPs and the genes in the region.
#' @references Hahne F. and Ivanek R. (2016). "Statistical Genomics: Methods and Protocols." In Mathe E and Davis S (eds.), chapter 
#' Visualizing Genomic Data Using Gviz and Bioconductor, pp. 335-351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi: 
#' 10.1007/978-1-4939-3578-9_16, \url{http://dx.doi.org/10.1007/978-1-4939-3578-9_16}. 
#' @export
plot_snp_module <- function(cones, k, genome = "hg19") {
  
  check_installed("GenomeInfoDb", "plot_snp_module")
  check_installed("GenomicRanges", "plot_snp_module")
  check_installed("Gviz", "plot_snp_module")
  check_installed("IRanges", "plot_snp_module")
  
  module <- subset(cones, module %in% k)
  
  by(module, module$chr, function(snps2plot) {

    snpRange <- GenomicRanges::GRanges(seqnames = paste0("chr", snps2plot$chr), 
                        ranges = IRanges::IRanges(start = snps2plot$pos, 
                                                  end = snps2plot$pos) )
    GenomeInfoDb::genome(snpRange) <- genome
    tsnp <- Gviz::AnnotationTrack(snpRange, name = "SNP", stacking ="dense")

    tideo <- Gviz::IdeogramTrack(genome = GenomeInfoDb::genome(snpRange), 
                           chromosome = names(GenomeInfoDb::genome(snpRange)))
    tseq <- Gviz::GenomeAxisTrack()
    
    tbio <- Gviz::BiomartGeneRegionTrack(genome = genome, 
                                   chromosome = names(GenomeInfoDb::genome(snpRange)),
                                   start = head(snps2plot$pos, n = 1), 
                                   end = tail(snps2plot$pos, n = 1), 
                                   name = "Ensembl")
    
    Gviz::plotTracks(list(tideo, tseq, tbio, tsnp), showId = TRUE)
  })
  
  return(TRUE)

}

#' Ideogram of SConES results.
#' 
#' @description Create a circular ideogram of the \code{SConES} results using the circlize package (Gu et al., 2014).
#' 
#' @param cones Output from \code{SConES}.
#' @param genome Abbreviations of the genome to use: hg19 for human (default), mm10 for mouse, etc. Argument to be passed to 
#' \code{\link{circos.initializeWithIdeogram}} \code{species}.
#' @return A circular ideogram, including the manhattan plot, and the interactions between the selected SNPs.
#' @references Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014). circlize Implements and enhances circular visualization in R. 
#' Bioinformatics (Oxford, England), 30(19), 2811-2. \url{https://doi.org/10.1093/bioinformatics/btu393}
#' @importFrom igraph %>%
#' @export
plot_ideogram <- function(cones, genome = "hg19") {
  
  check_installed("circlize", "plot_ideogram")
  
  circlize::circos.initializeWithIdeogram(species = genome)
  
  bed <- map2bed(cones)
  bed$c <- cones$c
  bed$selected <- cones$selected
  # order to put the selected snps in fron in the plot
  bed <- bed[with(bed, order(selected)),]
  
  circlize::circos.genomicTrackPlotRegion(bed, 
                                ylim = c(0, 1.1 * max(bed$c, na.rm = T)), 
                                panel.fun = function(region, value, ...) {
                                  # color according to selection/non-selection
                                  col = ifelse(value[[2]], "orange", "gray70")
                                  circlize::circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
                                }, track.height = 0.3)
  
  # create links
  selected <- subset(cones, selected)
  
  regions <- by(selected, selected[,c("chr","module")], function(k) {
    data.frame(chr = paste0("chr", unique(k$chr)),
               module = unique(k$module),
               start = min(k$pos),
               end = max(k$pos))
  }) %>% do.call(rbind, .)
  
  links <- merge(regions, regions, by = "module")
  
  region1 <- subset(links, chr.x != chr.y, select = c("chr.x","start.x","end.x"))
  region2 <- subset(links, chr.x != chr.y, select = c("chr.y","start.y","end.y"))
  
  circlize::circos.genomicLink(region1, region2, col = sample(1:5, nrow(region1), replace = TRUE))
  
  circlize::circos.clear()
  
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