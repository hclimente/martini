getGinSettings <- function(...){

  settings = list(...)

  # unsigned int
  if (! "nParameters" %in% settings)
    settings["nParameters"] = 10;

  # unsigned int
  if (! "folds" %in% settings)
    settings["folds"] = 10

  # bool, VectorXd, VectorXd
  if (! "autoParameters" %in% settings & "lambdas" %in% settings & "etas" %in% settings) {
    settings["autoParameters"] = TRUE
    settings["lambdas"] = rep(0, settings["nParameters"])
    settings["etas"] = rep(0, settings["nParameters"])
  }

  # unsigned int
  if (! "test_statistic" %in% settings)
    settings["test_statistic"] = 0

  # unsigned int
  if (! "gridsearch_depth" %in% settings)
    settings["gridsearch_depth"] = 1

  # unsigned int
  if (! "selection_criterion" %in% settings)
    settings["selection_criterion"] = 1

  # double
  if (! "seed" %in% settings)
    settings["seed"] = 0

  # double
  if (! "selection_ratio" %in% settings)
    settings["selection_ratio"] = 0.8

  # bool
  if (! "dump_intermediate_results" %in% settings)
    settings["dump_intermediate_results"] = TRUE

  # bool
  if (! "evaluateObjective" %in% settings)
    settings["evaluateObjective"] = FALSE

  return(settings);
}
