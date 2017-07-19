#' Ideogram of a SNP cluster.
#' 
#' @description Create a ideogram of a SNP cluster from \code{shake} results using the Gviz package (Hahne and Ivanek, 2016).
#' 
#' @param map Results from \code{shake}.
#' @param k Id of the cluster to plot.
#' @param genome Genome to use (by default hg19).
#' @return An ideogram per chromosome showing the selected SNPs and the genes in the region.
#' @references Hahne F. and Ivanek R. (2016). “Statistical Genomics: Methods and Protocols.” In Mathé E and Davis S (eds.), chapter Visualizing Genomic Data Using Gviz and Bioconductor, pp. 335-351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi: 10.1007/978-1-4939-3578-9_16, \url{http://dx.doi.org/10.1007/978-1-4939-3578-9_16}. 
#' @export
plot_snp_cluster <- function(map, k, genome = "hg19") {
  
  if (!requireNamespace("Gviz", quietly = TRUE))
    stop("Gviz needed for this function to work. Please install it.", call. = FALSE)
  
  cluster <- subset(map, cluster == k)
  
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
                                        name="Ensembl")
    
    plotTracks(list(itrack, gtrack, biomTrack, snpTrack), showId = TRUE)
  })
  
  return(TRUE)
}