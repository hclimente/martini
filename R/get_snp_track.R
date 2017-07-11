get_snp_track <- function(map, selected, name = "SNP"){

  if (!requireNamespace("Gviz", quietly = TRUE)) {
    stop("Gviz needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  s <- which(selected)
  f <- head(s, n=1) - 5
  l <- tail(s, n=1) + 5
  
  snpPlot <- map[f:l,]
  snps <- selected[f:l]
  
  snpRange <- GRanges(seqnames = paste0("chr", snpPlot$chr), 
                      ranges = IRanges(start = snpPlot$gpos, 
                                       end = snpPlot$gpos) )
  genome(snpRange) <- "hg19"
  snpTrack <- AnnotationTrack(snpRange, name = name, stacking ="dense")
  feature(snpTrack) <- ifelse(snps, "Selected", "Unselected")
  
  return(snpTrack)
}