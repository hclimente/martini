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
snp2gene <- function(gwas, organism = 9606, flank = 0) {
  
  # get map in appropriate format
  map <- gwas$map
  colnames(map) <- c("chr", "snp", "cm", "gpos", "allele1", "allele2")
  map$chr <- gsub("[Cc]hr", "", map$chr)
  
  # convert taxid to ensembl species name e.g. human databases are hsapiens_*
  urlTaxonomy <- "https://rest.ensembl.org/taxonomy"
  query <- paste0(urlTaxonomy, "/id/", organism, "?content-type=application/json")
  organism <- GET(query)
  organism <- content(organism, type ="application/json", encoding="UTF-8")
  organism <- unlist(strsplit(organism$name, " "))
  organism <- tolower(paste0(substr(organism[1], 1,1), organism[2]))
  
  # create mart from ENSEMBL
  # consider vertebrates and plants
  ensembl <- tryCatch({
    datasetName <- paste0(organism, "_gene_ensembl")
    biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = datasetName)
  }, error = function(e){
    tryCatch({
      datasetName <- paste0(organism, "_eg_gene")
      biomaRt::useMart("plants_mart", host="plants.ensembl.org", dataset = datasetName)  
    }, error = function(e) {
      stop(e)
    })
  })
  
  snp2gene <- by(map, map$chr, function(snps) {
    
    # get genes in the range defined by the snps
    grange <- paste(unique(snps$chr), min(snps$gpos) - flank, max(snps$gpos) + flank, sep = ":")
    genes <- biomaRt::getBM(attributes = c("ensembl_gene_id","external_gene_name","start_position","end_position"),
                   filters = c("chromosomal_region", "biotype"),
                   values = list(chromosomal_region=grange, biotype="protein_coding"), 
                   mart = ensembl)
    
    if(nrow(genes) == 0) {
      return()
    }
      
    # add a buffer before and after the gene
    genes$start_position <- genes$start_position - flank
    genes$start_position <- ifelse(genes$start_position < 0, 0, genes$start_position)
    genes$end_position <- genes$end_position + flank
    
    # convert to genomic ranges and check overlaps
    isnps <- with(snps, IRanges::IRanges(gpos, width=1, names=snp))
    igenes <- with(genes, IRanges::IRanges(start_position, end_position, names=ensembl_gene_id))
    olaps <- IRanges::findOverlaps(isnps, igenes)
    s2g <- cbind(snps[S4Vectors::queryHits(olaps),], genes[S4Vectors::subjectHits(olaps),])
    s2g$gene <- ifelse(s2g$external_gene_name == "" | is.na(s2g$external_gene_name), s2g$ensembl_gene_id, s2g$external_gene_name)
    subset(s2g, select = c("snp", "gene"))
    
  })
  snp2gene <- do.call("rbind", snp2gene)
  
  return(snp2gene)
  
}