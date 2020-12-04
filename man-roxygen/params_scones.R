#' @param score Association score to measure association between 
#' genotype and phenotype. Possible values: chi2 (default), glm.
#' @param criterion String with the function to measure the quality of a split. 
#' Possible values: consistency (default), bic, aic, aicc.
#' @param etas Numeric vector with the etas to explore in the grid search. If 
#' ommited, it's automatically created based on the association
#' scores.
#' @param lambdas Numeric vector with the lambdas to explore in the grid search.
#' If ommited, it's automatically created based on the association scores.
