#ifndef RSCONES2_TEST_ASSOCIATIONS
#define RSCONES2_TEST_ASSOCIATIONS

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "gwas/CScones.h"
#include "io/CSconesIO.h"
#include "io/CPlinkParser.h"
#include "globals.h"
#include "get_settings.h"

using namespace Rcpp;

// [[Rcpp::export]]
List run_shake(Eigen::MatrixXd X, Eigen::VectorXd Y, Eigen::SparseMatrix<double,Eigen::ColMajor> network, Rcpp::List userSettings) {

  CSconesSettings settings;
  settings = get_settings(userSettings);

  CScones scones;
  scones = CScones(Y, X, network, settings);

  cout << "Testing associations.\n";
  scones.test_associations();

  cout << "Collecting results.\n";
  VectorXd indicator = scones.getIndicatorVector();
  double objectiveScore = scones.getObjectiveScore();
  VectorXd scores = scones.getScoreStatistic();
  SparseMatrixXd W = scones.getW();

  return Rcpp::List::create(Rcpp::Named("indicator") = indicator,
                            Rcpp::Named("grid") = scones.getResultStack(),
                            Rcpp::Named("objective") = objectiveScore,
                            Rcpp::Named("scores") = scores,
                            Rcpp::Named("W") = W,
                            Rcpp::Named("lambda") = scones.getBestLambda(),
                            Rcpp::Named("eta") = scones.getBestEta());
}

#endif //RSCONES2_TEST_ASSOCIATIONS
