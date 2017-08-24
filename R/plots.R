#' Ideogram of a SNP cluster.
#' 
#' @description Create a ideogram of a SNP cluster from \code{shake} results using the Gviz package (Hahne and Ivanek, 2016).
#' 
#' @param cones Results from \code{shake}.
#' @param k Id of the cluster to plot.
#' @param genome Abbreviations of the genome to use: hg19 for human (default),  mm10 for mouse, etc. Argument to be passed to \code{\link{IdeogramTrack}}.
#' @return An ideogram per chromosome showing the selected SNPs and the genes in the region.
#' @references Hahne F. and Ivanek R. (2016). "Statistical Genomics: Methods and Protocols." In Mathe E and Davis S (eds.), chapter Visualizing Genomic Data Using Gviz and Bioconductor, pp. 335-351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi: 10.1007/978-1-4939-3578-9_16, \url{http://dx.doi.org/10.1007/978-1-4939-3578-9_16}. 
#' @importFrom GenomicRanges GRanges
#' @importFrom Gviz AnnotationTrack IdeogramTrack GenomeAxisTrack BiomartGeneRegionTrack plotTracks
#' @export
plot_snp_cluster <- function(cones, k, genome = "hg19") {
  
  cluster <- subset(cones, cluster == k)
  
  by(cluster, cluster$chr, function(snps2plot) {

    snpRange <- GRanges(seqnames = paste0("chr", snps2plot$chr), 
                        ranges = IRanges(start = snps2plot$pos, 
                                         end = snps2plot$pos) )
    genome(snpRange) <- genome
    snpTrack <- AnnotationTrack(snpRange, name = "SNP", stacking ="dense")

    itrack <- IdeogramTrack(genome = genome(snpRange), chromosome = names(genome(snpRange)))
    gtrack <- GenomeAxisTrack()
    
    biomTrack <- BiomartGeneRegionTrack(genome = genome, 
                                        chromosome = names(genome(snpRange)),
                                        start = head(snps2plot$pos, n = 1), 
                                        end = tail(snps2plot$pos, n = 1), 
                                        name = "Ensembl")
    
    plotTracks(list(itrack, gtrack, biomTrack, snpTrack), showId = TRUE)
  })
  
  return(TRUE)
}

#' Ideogram of shake results.
#' 
#' @description Create a circular ideogram of the \code{shake} results using the circlize package (Gu et al., 2014).
#' 
#' @param cones Output from \code{shake}.
#' @param net The SNP network provided to \code{shake}.
#' @param genome Abbreviations of the genome to use: hg19 for human (default), mm10 for mouse, etc. Argument to be passed to \code{\link{circos.initializeWithIdeogram}} \code{species}.
#' @return A circular ideogram, including the manhattan plot, and the interactions between the selected SNPs.
#' @references Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014). circlize Implements and enhances circular visualization in R. Bioinformatics (Oxford, England), 30(19), 2811-2. \url{https://doi.org/10.1093/bioinformatics/btu393}
#' @importFrom circlize circos.initializeWithIdeogram circos.genomicTrackPlotRegion circos.genomicPoints circos.genomicLink circos.clear
#' @importFrom igraph get.data.frame delete_vertices
#' @export
plot_ideogram <- function(cones, net, genome = "hg19"){
  
  circos.initializeWithIdeogram(species = genome)
  
  bed <- map2bed(cones)
  bed$c <- cones$c
  bed$selected <- cones$selected
  # order to put the selected snps in fron in the plot
  bed <- bed[with(bed, order(selected)),]
  
  circos.genomicTrackPlotRegion(bed, 
                                ylim = c(0, 1.1 * max(bed$c)), 
                                panel.fun = function(region, value, ...) {
                                  # color according to selection/non-selection
                                  col = ifelse(value[[2]], "orange", "gray70")
                                  circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
                                }, track.height = 0.3)
  
  # create links
  # remove unselected SNPs from the graph
  df <- delete_vertices(net, as.character(cones$snp[! cones$selected])) %>%
    get.data.frame %>%
    subset(select = c("from", "to"))
  
  region1 <- merge(df, cones, by.x = "from", by.y = "snp", sort = F) %>%
    subset(select = c("chr", "from", "cm", "pos")) %>%
    map2bed
  region2 <- merge(df, cones, by.x = "to", by.y = "snp", sort = F) %>%
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