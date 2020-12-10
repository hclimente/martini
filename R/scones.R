#' Find connected explanatory SNPs.
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being
#' connected in an underlying network. Select the hyperparameters by
#' cross-validation.
#' @template params_gwas
#' @template params_net
#' @template params_covars
#' @template params_score
#' @template params_criterion
#' @template params_scones
#' @template return_cones
#' @template reference_azencott
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' scones.cv(minigwas, gi)
#' scones.cv(minigwas, gi, score = "glm")
#' @export
scones.cv <- function(gwas, net, covars = data.frame(), score = "chi2", 
                      criterion = "consistency", etas = numeric(), 
                      lambdas = numeric()) {
  
  # set options
  opts <- parse_scones_settings(c = 1, score, criterion, etas, lambdas)
  c <- single_snp_association(gwas, covars, opts[['score']])
  opts <- parse_scones_settings(c = c, score, criterion, etas, lambdas)
  
  return(mincut.cv(gwas, net, covars, opts))
  
}

#' Compute Laplacian matrix
#' 
#' @template params_gwas
#' @template params_net
#' @return A Laplacian matrix.
#' @importFrom igraph simplify as_adj
#' @importFrom Matrix diag rowSums
#' @keywords internal
get_L <- function(gwas, net) {
  
  map <- sanitize_map(gwas)
  
  # remove redundant edges and self-edges in network and sort
  net <- simplify(net)
  L <- as_adj(net, type="both", sparse = TRUE, attr = "weight")
  L <- L[map[['snp']], map[['snp']]]
  L <- -L
  diag(L) <- rowSums(abs(L))
  
  return(L)
  
}

#' Run the cross-validated min-cut algorithm
#'
#' @template params_gwas
#' @template params_net
#' @template params_covars
#' @importFrom Matrix diag
#' @importFrom utils capture.output
#' @keywords internal
mincut.cv <- function(gwas, net, covars, opts) {
  
  L <- get_L(gwas, net)
  
  # prepare data: remove redundant edges and self-edges in network and sort
  gwas <- permute_snpMatrix(gwas)
  covars <- arrange_covars(gwas, covars) # TODO use PC as covariates
  
  cones <- sanitize_map(gwas)
  cones[['c']] <- single_snp_association(gwas, covars, opts[['score']])
  
  # grid search
  y <- gwas[['fam']][['affected']]
  
  K <- cut(seq(1, length(y)), breaks = 10, labels = FALSE)
  scores <- lapply(unique(K), function(k) {
    single_snp_association(gwas, covars, opts[['score']], 
                           samples = (K != k) )
  })
  
  grid_scores <- lapply(opts[['etas']], function(eta){
    lapply(opts[['lambdas']], function(lambda){
      folds <- lapply(scores, function(c) {
        c <- if (opts[['sigmod']]) c + lambda * diag(L) else c
        mincut_c(c, eta, lambda, L)
      })
      folds <- do.call(rbind, folds)
      score_fold(folds, opts[['criterion']], K, gwas, covars)
    }) %>% unlist
  })
  grid_scores <- do.call(rbind, grid_scores)
  dimnames(grid_scores) <- list(opts[['etas']], opts[['lambdas']])
 
  message('Grid of ', opts[['criterion']], ' scores (etas x lambdas):\n')
  message(paste0(capture.output(grid_scores), collapse = "\n") )
  
  best <- which(grid_scores == max(grid_scores), arr.ind = TRUE)
  best_eta <- opts[['etas']][best[1, 'row']]
  best_lambda <- opts[['lambdas']][best[1, 'col']]
  
  message("Selected parameters:\neta =",best_eta,"\nlambda =",best_lambda,"\n")
  
  c <- if (opts[['sigmod']]) cones[['c']] + best_lambda * diag(L) else cones[['c']]
  selected <- mincut_c(c, best_eta, best_lambda, L)
  cones[['selected']] <- as.logical(selected)
  
  cones <- get_snp_modules(cones, net)
  
  return(cones)
  
}

#' Find connected explanatory SNPs
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being
#' connected in an underlying network.
#' @template params_gwas
#' @template params_net
#' @param eta Value of the eta parameter.
#' @param lambda Value of the lambda parameter.
#' @template params_covars
#' @template params_score
#' @template return_cones
#' @inherit scones.cv return references
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' scones(minigwas, gi, 10, 1)
#' @export
scones <- function(gwas, net, eta, lambda, covars = data.frame(), score = 'chi2') {
  
  return(mincut(gwas, net, covars, eta, lambda, score, FALSE))
  
}

#' Run min-cut algorithm
#'
#' @template return_cones
#' @keywords internal
mincut <- function(gwas, net, covars, eta, lambda, score, sigmod) {
 
  L <- get_L(gwas, net)
   
  covars <- arrange_covars(gwas, covars) # TODO use PC as covariates
  
  cones <- sanitize_map(gwas)
  cones[['c']] <- single_snp_association(gwas, covars, score)
  c <- if (sigmod) cones[['c']] + lambda * diag(L) else cones[['c']]
     
  # run scores
  selected <- mincut_c(c, eta, lambda, L)
  cones[['selected']] <- as.logical(selected)
  
  cones <- get_snp_modules(cones, net)
  
  return(cones)
  
}

#' Calculate genotype-phhenotype associations 
#' 
#' @description Calculate the association between genotypes and a phenotype,
#' adjusting by covariates.
#' @template params_gwas
#' @template params_covars
#' @template params_score
#' @return A named vector with the association scores.
#' @importFrom snpStats single.snp.tests chi.squared snp.rhs.tests
#' @keywords internal
single_snp_association <- function(gwas, covars, score, 
                                   samples = rep(TRUE, nrow(gwas[['fam']])) ) {
  
  genotypes <- gwas[['genotypes']][samples, ]
  phenotypes <- gwas[['fam']][['affected']][samples]
  covars <- covars[samples, ]
    
  if (score == 'chi2') {
    
    tests <- single.snp.tests(phenotypes, snp.data = genotypes)
    c <- chi.squared(tests, df=1)
    
  } else if (score == 'glm') {
    
    if (ncol(covars) && nrow(covars) == nrow(genotypes)) {
      covars <- as.matrix(covars)
      tests <- snp.rhs.tests(phenotypes ~ covars, snp.data = genotypes)
    } else {
      tests <- snp.rhs.tests(phenotypes ~ 1, snp.data = genotypes)
    }
    
    c <- chi.squared(tests)
    names(c) <- names(tests)
  }
  
  c[is.na(c)] <- 0
  
  return(c)
}

#' Score the solutions of a k-fold
#' 
#' @description Take the k-solutions for a combination of hyperparameters, and 
#' assign a score to it (the larger, the better).
#' @param folds k-times-d matrix, where k is the number of folds, and d the 
#' number of SNPs.
#' @template params_criterion
#' @param K Numeric vector of length equal to the number of samples. The 
#' elements are integers from 1 to # folds. Indicates which samples belong to 
#' which fold.
#' @template params_gwas
#' @template params_covars
#' @importFrom methods as
#' @importFrom stats glm BIC AIC
#' @keywords internal
score_fold <- function(folds, criterion, K, gwas, covars) {
  
  score <- 0
  
  # penalize trivial solutions
  if (!sum(folds) | ( sum(folds)/length(folds) > .5) ){
    score <- -Inf
  } else if (criterion == 'consistency') {
    
    folds_diff <- unique(K)
    
    for (i in folds_diff) {
      for (j in folds_diff[folds_diff > i]) {
        C <- sum(folds[i,] & folds[j,])
        maxC <- sum(folds[i,] | folds[j,])
        score <- score + ifelse(maxC == 0, 0, C/maxC)
      }
    } 
  } else if (criterion %in% c('bic', 'aic', 'aicc')) {
    for (k in unique(K)) {
      
      phenotypes <- gwas[['fam']][['affected']][K != k]
      genotypes <- as(gwas[['genotypes']][K != k, ], 'numeric')
      genotypes <- as.data.frame(genotypes[ , folds[k,]])
      genotypes <- cbind(genotypes, covars[K != k, ])
      
      model <- glm(phenotypes ~ 1, data = genotypes)
      
      if (criterion == 'bic') {
        score <- score + BIC(model)
      } else if (criterion == 'aic') {
        score <- score + AIC(model)
      } else if (criterion == 'aicc') {
        k <- ncol(genotypes)
        n <- nrow(genotypes)
        score <- score + AIC(model) + (2*k^2 + 2*k)/(n - k - 1)
      }
    }
    # invert scores, as best models have the lowest scores
    score <- -score
  }
  
  return(score)
  
}

#' Return groups of interconnected SNPs.
#' 
#' @description Find modules composed by interconnected SNPs.
#' 
#' @param cones Results from \code{scones.cv}.
#' @template params_net
#' @return A list with the modules of selected SNPs.
#' @importFrom igraph induced_subgraph components
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' cones <- scones.cv(minigwas, gi)
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

#' Parse \code{scones.cv} settings
#' 
#' @description Creates a list composed by all \code{scones.cv} settings, with 
#' the values provided by the user, or the default ones if none is provided.
#' @template params_c
#' @template params_score
#' @template params_criterion
#' @template params_scones
#' @param sigmod Logical value stating if the settings are for SigMod.
#' @return A list of \code{scones.cv} settings.
#' @examples 
#' martini:::parse_scones_settings(etas = c(1,2,3), lambdas = c(4,5,6))
#' martini:::parse_scones_settings(c = c(1,10,100), score = "glm")
#' @keywords internal
parse_scones_settings <- function(c = numeric(), score = "chi2", 
                                  criterion = "consistency", etas = numeric(), 
                                  lambdas = numeric(), sigmod = FALSE) {
  
  settings <- list()
  
  valid_score <- c('glm', 'chi2')
  if (score %in% valid_score) {
    settings[['score']] <- score
  } else  {
    stop(paste("Error: invalid score", score))
  }
  
  valid_criterion <- c('consistency', 'bic', 'aic', 'aicc')
  if (criterion %in% valid_criterion) {
    settings[['criterion']] <- criterion
  } else {
    stop(paste("Error: invalid criterion", criterion))
  }
  
  logc <- log10(c[c != 0])
  if (length(etas) & is.numeric(etas)) {
    settings[['etas']] <- sort(etas)
  } else if (length(logc)) {
    minc <- min(logc)
    maxc <- max(logc)
    settings[['etas']] <- 10^seq(minc, maxc, length=10)
    settings[['etas']] <- signif(settings[['etas']], 3)
  } else {
    stop("Error: specify a valid etas or an association vector.")
  }
  
  if (length(lambdas) & is.numeric(lambdas)) {
    settings[['lambdas']] <- sort(lambdas)
  } else if (length(logc)) {
    minc <- min(logc)
    maxc <- max(logc)
    settings[['lambdas']] <- 10^seq(minc - 1, maxc + 1, length=10)
    settings[['lambdas']] <- signif(settings[['lambdas']], 3)
  } else {
    stop("Error: specify a valid lambdas or an association vector.")
  }
  
  if (length(sigmod) == 1 & is.logical(sigmod)) {
    settings[['sigmod']] <- sigmod
  } else {
    stop("Error: specify a logical sigmod.")
  }
  
  return(settings);
  
}