#' Simulate phenotype
#' 
#' @description Simulates a phenotype from a GWAS experiment and a specified set of causal SNPs.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snps Character vector with the SNP ids of the causal SNPs. Must match SNPs in gwas$map$snp.names.
#' @param h2 Heritability of the phenotype (between 0 and 1).
#' @param model String specifying the genetic model under the phenotype. Accepted values: "additive".
#' @param effectSize Numeric vector with the same lenght as the number of causal SNPs. It indicates the effect size of each of the SNPs; if absent, they are sampled fron a normal distribution.
#' @param qualitative Bool indicating if the phenotype is qualitative or not (quantitative).
#' @param ncases Integer specifying the number of cases to simulate in a qualitative phenotype. Required if qualitative = TRUE.
#' @param ncontrols Integer specifying the number of controls to simulate in a qualitative phenotype. Required if qualitative = TRUE.
#' @return A vector with the simulated phenotype for each sample.
#' @references Inspired from GCTA simulation tool: \url{http://cnsgenomics.com/software/gcta/Simu.html}.
#' @export
simulate_phenotype <- function(gwas, snps, h2, model = "additive", effectSize = rnorm(length(snps)), qualitative = FALSE, ncases, ncontrols){
  # TODO check correspondence with gcta implementation
  
  X <- as(gwas$genotypes, "numeric")
  
  if (any(! names(snps) %in% gwas$map$snp.names))
    stop(paste("The following causal SNPs are not in the SNP list:", setdiff(names(snps), gwas$map$snp.names)))
  if (h2 < 0 | h2 > 1)
    stop(paste0("h2 must be between 0 and 1. Current value is ", h2, "."))
  
  X <- X[, gwas$map$snp.names %in% names(snps)]
  
  # get effect sizes u and weights w
  u <- effectSize
  p <- (2 * colSums(X == 2) + colSums(X == 1)) / (2 * nrow(X))
  x <- 2 * (X == 2) + (X == 1)
  w = (x - 2 * p) / sqrt(2 * p * (1 - p))
  
  if (model == "additive")
    geno <- colSums(t(w) * u)
  else
    stop(paste0("Genetic model ", model, " not recognised."))
  
  residual.var <- var(geno) * (1 / h2 - 1)
  residual <- rnorm(length(geno), sd = sqrt(residual.var))
  
  trait <- geno + residual
  
  if (qualitative){
    if (! exists("ncases") )
      stop("Specify ncases if qualitative = TRUE.")
    else if (! exists("ncontrols") )
      stop("Specify ncontrols if qualitative = TRUE.")
    else if ( length(trait) < (ncases + ncontrols) )
      stop("Cases and controls requested exceed number of samples provided.")

    trait.sorted <- sort(trait, index.return = TRUE)
    cases <- head(trait.sorted$ix, n = ncases)
    controls <- tail(trait.sorted$ix, n = ncontrols)
    
    Y <- numeric(length(trait))
    Y <- NA
    Y[cases] <- 2
    Y[controls] <- 1
  } else {
    Y <- trait
  }
  
  return(Y)
}