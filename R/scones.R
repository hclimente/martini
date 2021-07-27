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
#' @template params_etas
#' @template params_lambdas
#' @template return_cones
#' @template params_family
#' @template params_link
#' @template reference_azencott
#' @template params_max_prop_snp
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' scones.cv(minigwas, gi)
#' scones.cv(minigwas, gi, score = "glm")
#' @export
scones.cv <- function(gwas, net, covars = data.frame(), 
                      score = c("chi2", "glm", "r2"), 
                      criterion = c("stability", "bic", "aic", "aicc", 
                                    "global_clustering", "local_clustering"), 
                      etas = numeric(), lambdas = numeric(),
                      family = c("binomial", "poisson", "gaussian", "gamma"), 
                      link = c("logit", "log", "identity", "inverse"),
                      max_prop_snp = 0.5) {
  
  # set options
  score <- match.arg(score)
  criterion <- match.arg(criterion)
  family <- match.arg(family)
  link <- match.arg(link)
  c <- snp_test(gwas, covars, score, family, link)
  grid <- get_grid(c = c, etas, lambdas)
  
  return(mincut.cv(gwas, net, covars, grid[['etas']], grid[['lambdas']], 
                   criterion, score, FALSE, family, link, max_prop_snp))
  
}

#' Run the cross-validated min-cut algorithm
#'
#' @template params_gwas
#' @template params_net
#' @template params_covars
#' @template params_family
#' @template params_link
#' @importFrom utils capture.output
#' @keywords internal
mincut.cv <- function(gwas, net, covars, etas, lambdas, criterion, score, 
                      sigmod, family, link, max_prop_snp) {
  
  # prepare data
  gwas <- permute_snpMatrix(gwas)
  A <- get_adjacency(gwas, net)
  
  # grid search
  K <- cut(seq(1, nrow(gwas[['fam']])), breaks = 10, labels = FALSE)
  folds <- lapply(unique(K), function(k) {
    
    gwas_k <- subset_snpMatrix(gwas, (K!=k))
    covars_k <- arrange_covars(gwas_k, covars)
    c_k <- snp_test(gwas_k, covars_k, score, family, link)
    
    lapply(lambdas, function(lambda) {
      c_k <- if (sigmod) c_k + lambda * rowSums(A) else c_k
      lapply(etas, function(eta) {
        selected_k <- mincut_c(c_k, eta, lambda, A)
        score_fold(gwas_k, covars_k, net, selected_k, criterion, max_prop_snp)
      })
    })
  })
  
  grid <- matrix(0, nrow = length(lambdas), ncol = length(etas), 
                 dimnames = list(lambdas, etas))
  
  for (i in seq(length(lambdas))) {
    for (j in seq(length(etas))){
      mat <- lapply(folds, function(x) x[[i]][[j]])
      mat <- do.call(rbind, mat)
      if (all(is.finite(mat))) {
        pearson_cor = cor(t(mat))
        grid[i,j] <- mean(pearson_cor[lower.tri(pearson_cor)], na.rm = TRUE)
      } else {
        grid[i, j] <- -Inf
      }
    }
  }
  
  best <- which(grid == max(grid), arr.ind = TRUE)
  best_lambda <- lambdas[best[1, 'row']]
  best_eta <- etas[best[1, 'col']]
  
  message('Grid of ', criterion, ' scores (lambdas \u00D7 etas):\n')
  message(paste0(capture.output(grid), collapse = "\n") )
  message("Selected parameters:\neta =", best_eta, "\nlambda =", best_lambda)
  
  cones <- mincut(gwas, net, covars, best_eta, best_lambda, score, sigmod, 
                  family, link)
  
  return(cones)
  
}

#' Find connected explanatory SNPs
#' 
#' @description Finds the SNPs maximally associated with a phenotype while being
#' connected in an underlying network.
#' @template params_gwas
#' @template params_net
#' @template params_net
#' @template params_eta
#' @template params_lambda
#' @template params_covars
#' @template params_score
#' @template params_family
#' @template params_link
#' @template return_cones
#' @inherit scones.cv return references
#' @examples
#' gi <- get_GI_network(minigwas, snpMapping = minisnpMapping, ppi = minippi)
#' scones(minigwas, gi, 10, 1)
#' @export
scones <- function(gwas, net, eta, lambda, covars = data.frame(), 
                   score = c("chi2", "glm", "r2"), 
                   family = c("binomial", "poisson", "gaussian", "gamma"), 
                   link = c("logit", "log", "identity", "inverse")) {
  
  score <- match.arg(score)
  family <- match.arg(family)
  link <- match.arg(link)
  
  return(mincut(gwas, net, covars, eta, lambda, score, FALSE, family, link))
  
}

#' Run min-cut algorithm
#'
#' @template return_cones
#' @keywords internal
mincut <- function(gwas, net, covars, eta, lambda, score, sigmod, family, link){
 
  A <- get_adjacency(gwas, net)
   
  covars <- arrange_covars(gwas, covars)
  
  map <- sanitize_map(gwas)
  c <- snp_test(gwas, covars, score, family, link)
  c <- if (sigmod) c + lambda * rowSums(A) else c
  selected <- mincut_c(c, eta, lambda, A)
  cones <- induced_subgraph(net, map[['snp']][selected])
  
  
  return(cones)
  
}

#' Calculate genotype-phhenotype associations 
#' 
#' @description Calculate the association between genotypes and a phenotype,
#' adjusting by covariates.
#' @template params_gwas
#' @template params_covars
#' @template params_score
#' @template params_family
#' @template params_link
#' @return A named vector with the association scores.
#' @importFrom snpStats single.snp.tests chi.squared snp.rhs.tests
#' @importFrom stats t.test
#' @keywords internal
snp_test <- function(gwas, covars, score, family, link) {
  
  genotypes <- gwas[['genotypes']]
  phenotypes <- gwas[['fam']][['affected']]
  
  if (score == 'chi2') {
    
    tests <- single.snp.tests(phenotypes, snp.data = genotypes)
    c <- chi.squared(tests, df=1)
    
  } else if (score == 'glm') {
    
    if (ncol(covars) && nrow(covars) == nrow(genotypes)) {
      covars <- as.matrix(covars)
      tests <- snp.rhs.tests(phenotypes ~ covars, snp.data = genotypes,
                             family = family, link = link)
    } else {
      tests <- snp.rhs.tests(phenotypes ~ 1, snp.data = genotypes,
                             family = family, link = link)
    }
    
    c <- chi.squared(tests)
    names(c) <- names(tests)
  } else if (score == 'r2') {
    
    c <- apply(genotypes, 2, cor, phenotypes, use="pairwise.complete.obs")
    c[is.na(c)] <- 0.01
    c <- c^2
    
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
#' @template params_max_prop_snp
#' (between 0 and 1). Larger solutions will be discarded.
#' @importFrom igraph induced_subgraph transitivity
#' @importFrom methods as
#' @importFrom stats glm BIC AIC
#' @keywords internal
score_fold <- function(gwas, covars, net, selected, criterion, max_prop_snp) {
  
  # score for a trivial solution
  score <- -Inf
  
  if (sum(selected) & (sum(selected) / length(selected) <= max_prop_snp) ){
    if (criterion == 'stability') {
      score <- selected
    } else if (criterion %in% c('bic', 'aic', 'aicc')) {
        
      phenotypes <- gwas[['fam']][['affected']]
      genotypes <- as(gwas[['genotypes']], 'numeric')
      genotypes <- as.data.frame(genotypes[, selected])
      if (ncol(covars)) {
        covars <- arrange_covars(gwas, covars)
        genotypes <- cbind(genotypes, covars)
      }
      
      model <- glm(phenotypes ~ ., data = genotypes)
      
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
      cones_subnet <- induced_subgraph(net, cones)
      
      if (criterion == 'local_clustering') {
        score <- transitivity(cones_subnet, type = 'local')
        score <- mean(score, na.rm = TRUE)
      } else if (criterion == 'global_clustering') {
        score <- transitivity(cones_subnet, type = 'global')
      }
    }
  }
  
  return(score)
  
}

#' Parse \code{scones.cv} settings
#' 
#' @description Creates a list composed by all \code{scones.cv} settings, with 
#' the values provided by the user, or the default ones if none is provided.
#' @template params_c
#' @template params_etas
#' @template params_lambdas
#' @return A list of \code{scones.cv} settings.
#' @importFrom stats quantile
#' @examples 
#' martini:::get_grid(etas = c(1,2,3), lambdas = c(4,5,6))
#' martini:::get_grid(c = c(1,10,100))
#' @keywords internal
get_grid <- function(c = numeric(), etas = numeric(), lambdas = numeric()) {
  
  grid <- list()
  
  if (length(c)) {
    # remove lowest c, which pull the grid towards irrelevant values
    maxc <- log10(max(c))
    ## add noise to avoid issues with repeated cs
    c <- c + abs(rnorm(length(c), sd = 1e-10))
    c <- c[c > quantile(c, 0.01)]
    minc <- log10(min(c))
  }
  
  if (length(etas) & is.numeric(etas)) {
    grid[['etas']] <- sort(etas)
  } else if (length(c)) {
    grid[['etas']] <- 10^seq(minc, maxc, length=10)
    grid[['etas']] <- signif(grid[['etas']], 3)
  } else {
    stop("Error: specify a valid etas or an association vector.")
  }
  
  if (length(lambdas) & is.numeric(lambdas)) {
    grid[['lambdas']] <- sort(lambdas)
  } else if (length(c)) {
    grid[['lambdas']] <- 10^seq(minc - 1, maxc + 1, length=10)
    grid[['lambdas']] <- signif(grid[['lambdas']], 3)
  } else {
    stop("Error: specify a valid lambdas or an association vector.")
  }
  
  return(grid)
  
}
