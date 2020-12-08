#' @template params_gwas
#' @template params_net
#' @template params_score
#' @template params_criterion
#' @param etas Numeric vector with the etas to explore in the grid search. If 
#' ommited, it's automatically created based on the association
#' scores.
#' @param lambdas Numeric vector with the lambdas to explore in the grid search.
#' If ommited, it's automatically created based on the association scores.
#' @template params_covars
