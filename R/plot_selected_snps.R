plot_selected_snps <- function(map){
  library(Gviz)
  
  s <- which(map$selected)
  f <- head(s, n=1) - 5
  l <- tail(s, n=1) + 5
  
  snpPlot <- map[f:l,]
  
  snpRange <- GRanges(seqnames = paste0("chr", snpPlot$chr), 
                      ranges = IRanges(start = snpPlot$gpos, 
                                       end = snpPlot$gpos) )
  genome(snpRange) <- "hg19"
  snpTrack <- AnnotationTrack(snpRange, name = "SNP", stacking ="dense")
  feature(snpTrack) <- ifelse(snpPlot$selected, "Selected", "Unselected")
  
  itrack <- IdeogramTrack(genome = genome(snpRange), chromosome = names(genome(snpRange)))
  gtrack <- GenomeAxisTrack()
  
  biomTrack <-BiomartGeneRegionTrack(genome = "hg19", 
                                     chromosome = names(genome(snpRange)),
                                     start = head(snpPlot$gpos, n = 1), 
                                     end = tail(snpPlot$gpos, n = 1), 
                                     name="Ensembl")
  
  plotTracks(list(itrack, gtrack, biomTrack, snpTrack), showId = TRUE,
             Selected = "red", Unselected = "grey20")
}