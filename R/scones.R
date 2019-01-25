#' Find connected explanatory SNPs.
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being
#'  connected in an underlying network (Azencott et al., 2013).
#' 
#' @param gwas A SnpMatrix object with the GWAS information.
#' @param net An igraph network that connects the SNPs.
#' @param encoding SNP encoding. Possible values: additive (default), resessive,
#' dominant, codominant.
#' @param covars A data frame with the covariates. It must contain a column 
#' 'sample' containing the sample IDs, and an additional columns for each 
#' covariate.
#' @param sigmod Boolean. If TRUE, use the Sigmod variant of SConES, meant to 
#' prioritize tightly connected clusters of SNPs
#' @param ... Extra arguments for \code{\link{parse_scones_settings}}.
#' @return A copy of the \code{SnpMatrix$map} \code{data.frame}, with the 
#' following additions:
#' \itemize{
#' \item{C: contains the univariate association score for every single SNP.}
#' \item{selected: logical vector indicating if the SNP was selected by SConES 
#' or not.}
#' \item{module: integer with the number of the module the SNP belongs to.}
#' }
#' @references Azencott, C. A., Grimm, D., Sugiyama, M., Kawahara, Y., & 
#' Borgwardt, K. M. (2013). Efficient network-guided multi-locus 
#' association mapping with graph cuts. Bioinformatics, 29(13), 171-179. 
#' \url{https://doi.org/10.1093/bioinformatics/btt238}
#' @importFrom igraph simplify as_adj
#' @importFrom Matrix diag rowSums
#' @importFrom methods as
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' search_cones(minigwas, gi)
#' search_cones(minigwas, gi, encoding = "recessive")
#' search_cones(minigwas, gi, associationScore = "skat")
#' @export
search_cones <- function(gwas, net, encoding = "additive", 
                         covars = data.frame(), sigmod = FALSE, ...) {

  cones <- gwas[["map"]]
  colnames(cones) <- c("chr","snp","cm","pos","allele.1", "allele.2")
  cones[['c']] <- single_snp_association(gwas[['genotypes']], 
                                         gwas[['fam']][['affected']], 
                                         covariates, 
                                         settings[['associationScore']])
  
  settings <- parse_scones_settings(c = cones[['c']], 
                                    modelScore = "consistency", ...)
  
  # prepare network
  ## remove redundant edges and self-edges
  net <- simplify(net)
  W <- as_adj(net, type="both", sparse = TRUE, attr = "weight")
  diag(W) <- ifelse(sigmod, -rowSums(W), 0)
  
  ## order according to order in map
  W <- W[cones[['snp']], cones[['snp']]]
  
  # K-fold association
  associationScores <- list()
  ids <- gwas[['fam']][['affected']]
  
  K <- 10
  K <- cut(seq(1, length(ids)), breaks=K,labels=FALSE)
  for (k in unique(K)) {
      
      genotypes <- gwas[['genotypes']][K != k, ]
      phenotypes <- gwas[['fam']][['affected']][K != k]
      covariates <- covars[K != k, ]
      
      associationScores[[k]] <- single_snp_association(genotypes, phenotypes, 
                                    covariates, settings[['associationScore']])

  }
  
  best_eta <- best_lambda <- best_score <- 0
  
  for (eta in settings$etas) {
      for (lambda in settings$lambdas) {
          
          # TODO if solution is empty, break
          
          folds <- lapply(associationScores, run_scones, eta, lambda, W)
          folds <- do.call(rbind, folds)
          
          if (score_fold(folds, settings[['modelScore']]) > best_score) {
              best_eta <- eta
              best_lambda <- lambda
          }
      }
  }

  cat("eta =", best_eta, "\nlambda =", best_lambda, "\n")
  
  selected <- run_scones(cones[['c']], best_eta, best_lambda, W)
  cones[['selected']] <- as.logical(selected)
  
  cones <- get_snp_modules(cones, net)
  
  return(cones)
  
}

#' @importFrom snpStats single.snp.tests chi.squared
single_snp_association <- function(genotypes, phenotypes, 
                                   covars, associationScore) {
    
    if (TRUE) {
        
        tests <- single.snp.tests(phenotypes, snp.data = genotypes)
        c <- chi.squared(tests, df=1)
        
    } else {
        X <- as(gwas[['genotypes']], "numeric")
        X <- encode_gwas(X, encoding)
        Y <- gwas[['fam']][['affected']]
        
    }
    
    c[is.na(c)] <- 0
    
    return(c)
}

score_fold <- function(folds, model_score) {
    
    k <- nrow(folds)
    score <- 0
    
    for (i in seq(k - 1)) {
        for (j in seq(i + 1, k)) {
            
            if (model_score == 'consistency') {
                C <- sum(folds[i,] * folds[j,])
                maxC <- max(sum(folds[i,]), sum(folds[j,]))
                score <- score + ifelse(maxC == 0, 0, C/maxC)
            }
        }
    }
    
    return(score)
    
}

#' Return groups of interconnected SNPs.
#' 
#' @description Find modules composed by interconnected SNPs.
#' 
#' @param cones Results from \code{evo}.
#' @param net The same SNP network provided to \code{evo}.
#' @return A list with the modules of selected SNPs.
#' @importFrom igraph induced_subgraph components
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' cones <- search_cones(minigwas, gi)
#' martini:::get_snp_modules(cones, gi)
#' @keywords internal
get_snp_modules <- function(cones, net) {
  
  selected <- subset(cones, selected)
  subnet <- induced_subgraph(net, as.character(selected[,'snp']))
  
  modules <- components(subnet)
  modules <- as.data.frame(modules['membership'])
  colnames(modules) <- "module"
  modules['snp'] <- rownames(modules)
  
  modules <- merge(cones, modules, all.x = TRUE)
  cones <- modules[match(cones[,'snp'], modules[,'snp']),]

  return(cones)
  
}

#' Parse search_cones settings.
#' 
#' @description Creates a list composed by all \code{evo} settings, with the 
#' values provided by the user, or the default ones if none is provided.
#' @param associationScore Association score to measure association between 
#' genotype and phenotype. Possible values: chi2 (default), skat, trend.
#' @param modelScore Model selection criterion Possible values: consistency, 
#' bic (default), aic, aicc, mbic.
#' @param etas Numeric vector with the etas to explore in the grid search. If 
#' ommited, it's automatically created based on the association
#' scores.
#' @param lambdas Numeric vector with the lambdas to explore in the grid search.
#' If ommited, it's automatically created based on the association scores.
#' @param debug Display additional information. Possible values: TRUE, FALSE
#' (default).
#' @param c Numeric vector with the association scores of the SNPs. Specify it 
#' to automatically an appropriate range of etas and lambas.
#' @return A list of \code{evo} settings.
#' @examples 
#' martini:::parse_scones_settings(etas = c(1,2,3), lambdas = c(4,5,6))
#' martini:::parse_scones_settings(c = c(1,10,100), associationScore = "skat")
#' @keywords internal
parse_scones_settings <- function(associationScore = "chi2", modelScore = "bic", 
                             etas = numeric(), lambdas = numeric(), 
                             debug = FALSE, c = numeric()){
  
  settings <- list()
  
  # unsigned int
  valid_associationScore <- c('skat', 'chi2')
  if (associationScore %in% valid_associationScore) {
    settings[['associationScore']] <- associationScore
  } else  {
    stop(paste("Error: invalid associationScore", associationScore))
  }
  
  # unsigned int
  valid_modelScore <- c('consistency', 'bic', 'aic', 'aicc', 'mbic')
  if (modelScore %in% valid_modelScore) {
    settings[['modelScore']] <- modelScore
  } else {
    stop(paste("Error: invalid modelScore", modelScore))
  }
  
  # bool
  if (! is.logical(debug)) {
    stop("Error: debug must be logical.")
  } else {
    settings[['debug']] <- debug;
  }
  
  # VectorXd
  logc <- log10(c[c != 0])
  minc <- min(logc)
  maxc <- max(logc)
  if (length(etas) & is.numeric(etas)) {
    settings[['etas']] <- sort(etas, decreasing = TRUE)
  } else if (length(logc)) {
    settings[['etas']] <- seq(maxc, minc, length=10)
  } else {
      stop("Error: specify a valid etas or an association vector.")
  }
  
  # VectorXd
  if (length(lambdas) & is.numeric(lambdas)) {
      settings[['lambdas']] <- sort(lambdas, decreasing = TRUE)
  } else if (length(logc)) {
      settings[['lambdas']] <- seq(maxc, minc, length=10)
  } else {
      stop("Error: specify a valid lambdas or an association vector.")
  }
  
  return(settings);
  
}

#' Get evo settings.
#' 
#' @description Creates a list composed by all \code{evo} settings, with the 
#' values provided by the user, or the default ones if none is provided.
#' @param associationScore Association score to measure association between 
#' genotype and phenotype. Possible values: chi2 (default), skat, trend.
#' @param modelScore Model selection criterion Possible values: consistency, 
#' bic (default), aic, aicc, mbic.
#' @param etas Numeric vector with the etas to explore in the grid search. If 
#' ommited, it's automatically created based on the association
#' scores.
#' @param lambdas Numeric vector with the lambdas to explore in the grid search.
#' If ommited, it's automatically created based on the association scores.
#' @param debug Display additional information. Possible values: TRUE, FALSE
#' (default).
#' @return A list of \code{evo} settings.
#' @examples 
#' martini:::get_evo_settings()
#' martini:::get_evo_settings(associationScore = "skat")
#' @keywords internal
get_evo_settings <- function(associationScore = "chi2", modelScore = "bic", 
                             etas = numeric(), lambdas = numeric(), 
                             debug = FALSE){
    
    settings <- list()
    
    # unsigned int
    settings[['associationScore']] <- switch(associationScore, skat = 0, chi2 = 1)
    if (length(settings[['associationScore']]) == 0) {
        stop(paste("Error: invalid associationScore", associationScore))
    }
    
    # unsigned int
    settings[['modelScore']] <- switch(modelScore, consistency = 0, bic = 1, 
                                       aic = 2, aicc = 3, mbic = 4)
    if (length(settings[['modelScore']]) == 0) {
        stop(paste("Error: invalid modelScore", modelScore))
    }
    
    # bool
    if (! is.logical(debug)) {
        stop("Error: debug must be logical.")
    } else {
        settings[['debug']] <- debug;
    }
    
    # VectorXd
    if (!is.numeric(etas)){
        stop("Error: etas must be numeric")
    } else {
        settings[['etas']] <- etas
    }
    
    # VectorXd
    if (!is.numeric(lambdas)){
        stop("Error: lambdas must be numeric")
    } else {
        settings[['lambdas']] <- lambdas
    }
    
    return(settings);
}