#' @return A copy of the \code{SnpMatrix$map} \code{data.frame}, with the 
#' following additions:
#' \itemize{
#' \item{c: contains the univariate association score for every single SNP.}
#' \item{selected: logical vector indicating if the SNP was selected by SConES 
#' or not.}
#' \item{module: integer with the number of the module the SNP belongs to.}
#' }
