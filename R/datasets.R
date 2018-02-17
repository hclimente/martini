#' Description of the minigwas dataset.
#'
#' @name minigwas
#' @docType data
#' @description Small GWAS example.
#' @format A list with 3 items:
#' \describe{
#'   \item{genotypes}{Genotype and phenotype information.}
#'   \item{fam}{Simulated network.}
#'   \item{map}{Result of runing \code{find_cones} with gwas and net.}
#' }
#' @examples 
#' data(minigwas)
#' 
#' # access different elements
#' minigwas[["genotypes"]]
#' minigwas[["map"]]
#' minigwas[["fam"]]
NULL

#' Genes for the minigwas dataset.
#'
#' @name minisnpMapping
#' @docType data
#' @description data.frame that maps SNPs from minigwas to their gene.
#' @examples
#' data(minisnpMapping)
#' 
#' head(minisnpMapping)
NULL

#' PPIs for the minigwas dataset.
#'
#' @name minippi
#' @docType data
#' @description data.frame describing pairs of proteins that interact for 
#' minigwas.
#' @examples
#' data(minippi)
#' 
#' head(minippi)
NULL