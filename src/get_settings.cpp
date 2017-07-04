#include "getSettings.h"

CSconesSettings getSettings(List userSettings) {

  CSconesSettings settings;
  settings = CSconesSettings();

  // unsigned int
  if (userSettings.containsElementNamed("folds"))
    settings.folds = userSettings["folds"];

  // bool, VectorXd, VectorXd
  if (userSettings.containsElementNamed("autoParameters") &
      userSettings.containsElementNamed("lambdas") &
      userSettings.containsElementNamed("etas")){
    settings.autoParameters = userSettings["autoParameters"];
    settings.lambdas = userSettings["lambdas"];
    settings.etas = userSettings["etas"];
  }

  // unsigned int
  if (userSettings.containsElementNamed("test_statistic"))
    settings.test_statistic = userSettings["test_statistic"];

  // unsigned int
  if (userSettings.containsElementNamed("nParameters"))
    settings.nParameters = userSettings["nParameters"];

  // unsigned int
  if (userSettings.containsElementNamed("gridsearch_depth"))
    settings.gridsearch_depth = userSettings["gridsearch_depth"];

  // unsigned int
  if (userSettings.containsElementNamed("selection_criterion"))
    settings.selection_criterion = userSettings["selection_criterion"];

  // double
  if (userSettings.containsElementNamed("seed"))
    settings.seed = userSettings["seed"];

  // double
  if (userSettings.containsElementNamed("selection_ratio"))
    settings.selection_ratio = userSettings["selection_ratio"];

  // bool
  if (userSettings.containsElementNamed("dump_intermediate_results"))
    settings.dump_intermediate_results = userSettings["dump_intermediate_results"];

  // bool
  if (userSettings.containsElementNamed("evaluateObjective"))
    settings.evaluateObjective = userSettings["evaluateObjective"];

  return settings;

}
