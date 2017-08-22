#' Get shake settings.
#' 
#' @description Creates a list composed by all \code{shake} settings, with the values provided by the user, or the default ones if none is provided.
#' 
#' @param ... Any \code{shake} option.
#' @return A list of \code{shake} settings.
get_evo_settings <- function(...){

  settings <- list(...)
  
  # unsigned int
  if (! "nParameters" %in% names(settings))
    settings[["nParameters"]] = 10;
  
  # unsigned int
  if (! "encoding" %in% names(settings))
    settings[["encoding"]] = 0;

  # unsigned int
  if (! "folds" %in% names(settings))
    settings[["folds"]] = 10

  # bool, VectorXd, VectorXd
  if (! ("autoParameters" %in% names(settings) & "lambdas" %in% names(settings) & "etas" %in% names(settings)) ) {
    settings[["autoParameters"]] = TRUE
    settings[["lambdas"]] = numeric()
    settings[["lambdas"]] = rep(0, settings[["nParameters"]])
    settings[["etas"]] = numeric()
    settings[["etas"]] = rep(0, settings[["nParameters"]])
  }

  # unsigned int
  if (! "associationScore" %in% names(settings))
    settings[["associationScore"]] = 0

  # unsigned int
  if (! "gridsearch_depth" %in% names(settings))
    settings[["gridsearch_depth"]] = 1

  # unsigned int
  if (! "modelScore" %in% names(settings))
    settings[["modelScore"]] = 1

  # double
  if (! "seed" %in% names(settings))
    settings[["seed"]] = 0

  # double
  if (! "selection_ratio" %in% names(settings))
    settings[["selection_ratio"]] = 0.8

  # bool
  if (! "dump_intermediate_results" %in% names(settings))
    settings[["dump_intermediate_results"]] = TRUE

  # bool
  if (! "evaluateObjective" %in% names(settings))
    settings[["evaluateObjective"]] = FALSE

  return(settings);
}
