#' Simulate causal SNPs
#' 
#' @description Selects randomly interconnected genes as causal, then selects a
#' proportion of them as causal.
#' 
#' @param net An igraph gene-interaction (GI) network that connects the SNPs.
#' @param ngenes Number of causal genes.
#' @param pcausal Number between 0 and 1, proportion of the SNPs in causal genes
#' that are causal themselves.
#' @return A vector with the ids of the simulated causal SNPs.
#' @importFrom igraph components V neighbors %>%
#' @importFrom stats na.omit
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' simulate_causal_snps(gi, ngenes=2)
#' @export
simulate_causal_snps <- function(net, ngenes = 20, pcausal = 1) {
  
  # SNPs with only one gene
  subnet <- subnet(net, "nGenes", 1)
  # genes with more than 1 SNP
  genes <- names(which(table(V(subnet)$gene) > 5))
  
  repeat {
    g <- sample(genes, 1)
    seed <- subvert(subnet, "gene", g)[1]
    
    causalGenes <- neighbors(subnet,seed)$gene %>% na.omit %>% unique
    causalGenes <- c(causalGenes, g)
    
    causalNetwork <- lapply(causalGenes, function(s) {
      subvert(subnet, "gene", s)$gene
    }) %>% unlist %>% table
    
    causalGenes <- causalNetwork[causalNetwork > ngenes/4] %>% names %>% 
                      intersect(genes) %>% head(ngenes)
    
    if ( length(causalGenes) >= ngenes ) {
      
      causal <- subvert(subnet, "gene", causalGenes)
      # sample a proportion pcausal of the snps in the causal genes
      causal <- sample(causal, length(causal) * pcausal)
      
      # check that we have at most two subnetworks
      if (components(subnet(net, "name", names(causal)))$no <= 2) {
        break
      }
    }
  }
  
  # get the nodes from the original net
  causal <- subvert(net, "name", names(causal))
  
  return(causal)
}

#' Simulate phenotype
#' 
#' @description Simulates a phenotype from a GWAS experiment and a specified set
#' of causal SNPs. If the data is qualitative, only controls are used.
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param snps Character vector with the SNP ids of the causal SNPs. Must match
#' SNPs in gwas[["map"]][["snp.names"]].
#' @param h2 Heritability of the phenotype (between 0 and 1).
#' @param model String specifying the genetic model under the phenotype.
#' Accepted values: "additive".
#' @param effectSize Numeric vector with the same lenght as the number of causal
#' SNPs. It indicates the effect size of each of the SNPs; if absent, they are
#' sampled fron a normal distribution.
#' @param qualitative Bool indicating if the phenotype is qualitative or not
#' (quantitative).
#' @param ncases Integer specifying the number of cases to simulate in a
#' qualitative phenotype. Required if qualitative = TRUE.
#' @param ncontrols Integer specifying the number of controls to simulate in a
#' qualitative phenotype. Required if qualitative = TRUE.
#' @param prevalence Value between 0 and 1 specifying the population prevalence
#' of the disease. Note that ncases cannot be greater than 
#' prevalence * number of samples. Required if qualitative = TRUE.
#' @return A copy of the GWAS experiment with the new phenotypes in
#' \code{gwas[["fam"]][["affected"]]}.
#' @references Inspired from GCTA simulation tool:
#' \url{http://cnsgenomics.com/software/gcta/Simu.html}.
#' @importClassesFrom snpStats SnpMatrix
#' @importFrom methods as
#' @importFrom stats rnorm var
#' @importFrom utils head tail
#' @examples 
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' causal <- simulate_causal_snps(gi, ngenes = 2)
#' simulate_phenotype(minigwas, causal, h2 = 1)
#' @export
simulate_phenotype <- function(gwas, snps, h2, model = "additive", 
                               effectSize = rnorm(length(snps)), 
                               qualitative = FALSE, 
                               ncases, ncontrols, prevalence){

  # select only controls
  binary <- (unique(gwas[["fam"]][["affected"]]) %>% length) == 2
  if (binary) {
    gwas[["genotypes"]] <- gwas[["genotypes"]][gwas[["fam"]][["affected"]]==1, ]
    gwas[["fam"]] <- gwas[["fam"]][gwas[["fam"]][["affected"]] == 1, ]
  }
  
  X <- as(gwas[["genotypes"]], "numeric")
  
  if (any(! names(snps) %in% gwas[["map"]][["snp.names"]])) {
    stop(paste("The following causal SNPs are not in the SNP list:", 
         paste(setdiff(names(snps), gwas[["map"]][["snp.names"]]), 
               collapse = ",")))
  }
  
  if (h2 < 0 | h2 > 1) {
    stop(paste0("h2 must be between 0 and 1. Current value is ", h2, "."))
  }
  
  X <- X[, names(snps)]
  G <- calculateG(effectSize, X, model)
  E <- calculateE(G, h2)
  Y <- G + E
  
  if (qualitative){
    if ( length(Y) < (ncases + ncontrols) ) {
      stop("Requested number of cases and controls too high (> # samples).")
    } else if ( prevalence * length(Y) < ncases ) {
      stop("Requested number of cases too high (> # samples * prevalence).")
    }
    
    Y.sorted <- sort(Y, index.return = TRUE)
    cases <- head(Y.sorted[['ix']], n = ncases)
    controls <- tail(Y.sorted[['ix']], n = ncontrols)

    Y[cases] <- 2
    Y[controls] <- 1
    Y[-c(cases,controls)] <- NA
    
  }
  
  gwas[["fam"]][["affected"]] <- Y
  
  return(gwas)
}

#' Calculate the genetic component of the phenotype
#' 
#' @description Calculates the genetic component of the phenotype from a
#' genotype.
#' 
#' @param effectSize A vector with the effect size of each SNP.
#' @param X Genotypes in a numeric matrix, where each row is a sample and each
#' column a SNP.
#' @param model Genetic model to assume.
#' @return A vector with the genetic component of each sample.
#' @keywords internal
calculateG <- function(effectSize, X, model) {
  
  X <- t(X)
  
  # calculate weights w
  p <- (2 * rowSums(X == 0) + rowSums(X == 1)) / (2 * ncol(X))
  x <- 2 * (X == 0) + (X == 1)
  w <- (x - 2 * p) / sqrt(2 * p * (1 - p))
  
  if (model == "additive") {
    G <- colSums(w * effectSize)
  } else {
    stop(paste0("Genetic model ", model, " not recognised."))
  }
  
  return(G)
  
}

#' Calculate the environmental component of the phenotype
#' 
#' @description Calculates the environmental component of the phenotype using
#' the variance in the genetic component.
#' 
#' @param G The genetic component of the phenotype.
#' @param h2 The heritability.
#' @return A vector with the environmental component of each sample.
#' @keywords internal
calculateE <- function(G, h2) {
  
  residual.var <- var(G) * (1/h2 - 1)
  residual <- rnorm(length(G), sd = sqrt(residual.var))
  
  return(residual)
  
}