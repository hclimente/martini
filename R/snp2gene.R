#' Map SNPs to genes.
#' 
#' @description Maps SNPs from a GWAS experiment to genes.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param organism Organism: hsapiens represents human, athaliana for Arabidopsis thaliana, etc.
#' @param flank A number with the flanking regions around genes to be considered part of the gene 
#' i.e. SNPs mapped to them will be considered mapped to the gene.
#' @return A dataframe with two columns: one for the SNP and another for the gene it has been 
#' mapped to.
#' @export
snp2gene <- function(gwas, organism = "hsapiens", flank = 0) {
  
  if (!requireNamespace("biomaRt", quietly = TRUE))
    stop("biomaRt needed for this function to work. Please install it.", call. = FALSE)
  if (!requireNamespace("IRanges", quietly = TRUE))
    stop("IRanges needed for this function to work. Please install it.", call. = FALSE)
  
  # get map in appropriate format
  map <- gwas$map
  colnames(map) <- c("chr", "snp", "cm", "gpos", "allele1", "allele2")
  map$chr <- gsub("[Cc]hr", "", map$chr)
  
  # create mart from ENSEMBL
  # consider vertebrates and plants
  ensembl <- tryCatch({
    datasetName <- paste0(organism, "_gene_ensembl")
    useMart("ENSEMBL_MART_ENSEMBL", dataset = datasetName)
  }, error = function(e){
    if (e$call == "useDataset(mart = mart, dataset = dataset, verbose = verbose)") {
      datasetName <- paste0(organism, "_eg_gene")
      useMart("plants_mart", host="plants.ensembl.org", dataset = datasetName)  
    } else {
      stop(e)
    }
  })
  
  by(map, map$chr, function(x) {
    grange <- paste(unique(x$chr), min(x$gpos), max(x$gpos), sep = ":")
    getBM(attributes = c("external_gene_name"),
          filters = c("chromosomal_region"),
          values = grange, 
          mart = ensembl)
  })
  
  snp2gene <- by(map, map$chr, function(snps) {
    
    # get genes in the range defined by the snps
    grange <- paste(unique(snps$chr), min(snps$gpos), max(snps$gpos), sep = ":")
    genes <- getBM(attributes = c("ensembl_gene_id","start_position", "end_position"),
                   filters = "chromosomal_region",
                   values = grange, 
                   mart = ensembl)
    
    # add a buffer before and after the gene
    chrEnd <- max(genes$end_position)
    genes$start_position <- genes$start_position - flank
    genes$start_position <- ifelse(genes$start_position < 0, 0, genes$start_position)
    genes$end_position = genes$end_position + flank
    genes$end_position <- ifelse(genes$end_position > chrEnd, chrEnd, genes$end_position)
    
    # convert to genomic ranges and check overlaps
    isnps <- with(snps, IRanges(gpos, width=1, names=snp))
    igenes <- with(genes, IRanges(start_position, end_position, names=ensembl_gene_id))
    olaps <- findOverlaps(isnps, igenes)
    s2g <- cbind(snps[queryHits(olaps),], genes[subjectHits(olaps),])
    subset(s2g, select = c("snp", "ensembl_gene_id"))
    
  })
  snp2gene <- do.call("rbind", snp2gene)
  
  return(snp2gene)
  
}