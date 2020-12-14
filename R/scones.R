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
  opts[['gwas']] <- gwas
  opts[['net']] <- net
  opts[['covars']] <- covars
  
  return(do.call(mincut.cv, opts))
  
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
mincut.cv <- function(gwas, net, covars, lambdas, etas, criterion, score, sigmod) {
  
  # prepare data
  gwas <- permute_snpMatrix(gwas)
  covars <- arrange_covars(gwas, covars) # TODO use PC as covariates
  L <- get_L(gwas, net)
  
  # grid search
  K <- cut(seq(1, nrow(gwas[['fam']])), breaks = 10, labels = FALSE)
  folds <- lapply(unique(K), function(k) {
    
    gwas_k <- subset_snpMatrix(gwas, (K!=k))
    covars_k <- covars[(K!=k), ]
    c_k <- single_snp_association(gwas_k, covars_k, score)
    
    lapply(lambdas, function(lambda) {
      c_k <- if (sigmod) c_k + lambda * diag(L) else c_k
      lapply(etas, function(eta) {
        selected_k <- mincut_c(c_k, eta, lambda, L)
        score_fold(gwas_k, covars_k, net, selected_k, criterion)
      })
    })
  })
  
  grid <- matrix(0, nrow = length(lambdas), ncol = length(etas), 
                 dimnames = list(lambdas, etas))
  
  for (i in seq(length(lambdas))) {
    for (j in seq(length(etas))){
      mat <- lapply(folds, function(x) x[[i]][[j]])
      mat <- do.call(rbind, mat)
      grid[i,j] <- sum(colSums(mat)) / (10 * sum(colSums(mat) != 0))
    }
  }
  
  best <- which(grid == max(grid), arr.ind = TRUE)
  best_lambda <- lambdas[best[1, 'row']]
  best_eta <- etas[best[1, 'col']]
  
  message('Grid of ', criterion, ' scores (lambdas×etas):\n')
  message(paste0(capture.output(grid), collapse = "\n") )
  message("Selected parameters:\neta =", best_eta, "\nlambda =", best_lambda)
  
  cones <- sanitize_map(gwas)
  cones[['c']] <- single_snp_association(gwas, covars, score)
  c <- if (sigmod) cones[['c']] + best_lambda * diag(L) else cones[['c']]
  cones[['selected']] <- mincut_c(c, best_eta, best_lambda, L)
  
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
  cones[['selected']] <- mincut_c(c, eta, lambda, L)
  
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
  
  genotypes <- gwas[['genotypes']]
  phenotypes <- gwas[['fam']][['affected']]
  
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
#' @template params_gwas
#' @template params_net
#' @template params_covars
#' @template params_criterion
#' @param max_soluction Maximum fraction of the SNPs involved in the solution
#' (between 0 and 1). Larger solutions will be discarded.
#' @importFrom igraph induced.subgraph transitivity
#' @importFrom methods as
#' @importFrom stats glm BIC AIC
#' @keywords internal
score_fold <- function(gwas, covars, net, selected, criterion, max_solution = .5) {
  
  # score for a trivial solution
  score <- -Inf
  
  if (sum(selected) & (sum(selected)/length(selected) <= max_solution) ){
    if (criterion == 'consistency') {
      score <- selected
    } else if (criterion %in% c('bic', 'aic', 'aicc')) {
        
      phenotypes <- gwas[['fam']][['affected']]
      genotypes <- as(gwas[['genotypes']], 'numeric')
      genotypes <- as.data.frame(genotypes[, selected])
      if (ncol(covars)) genotypes <- cbind(genotypes, covars)
      
      model <- glm(phenotypes ~ 1, data = genotypes)
      
      if (criterion == 'bic') {
        score <- BIC(model)
      } else if (criterion == 'aic') {
        score <- AIC(model)
      } else if (criterion == 'aicc') {
        k <- ncol(genotypes)
        n <- nrow(genotypes)
        score <- AIC(model) + (2*k^2 + 2*k)/(n - k - 1)
      }
      # invert scores, as best models have the lowest scores
      score <- -score
    } else if (criterion %in% c('local_clustering', 'global_clustering')) {
      cones <- sanitize_map(gwas)
      cones <- cones[selected, 'snp']
      cones_subnet <- induced.subgraph(net, cones)
      
      if (criterion == 'local_clustering') {
        score <- transitivity(cones_subnet, type = 'local')
        score <- mean(score, na.rm = T)
      } else if (criterion == 'global_clustering') {
        score <- transitivity(cones_subnet, type = 'global')
      }
    }
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
  
  valid_criterion <- c('consistency', 'bic', 'aic', 'aicc', 'global_clustering',
                       'local_clustering')
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
