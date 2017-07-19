#ifndef SHAKE_TEST_ASSOCIATIONS
#define SHAKE_TEST_ASSOCIATIONS

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/gwas/CScones.h"
#include "gin/io/CSconesIO.h"
#include "gin/io/CPlinkParser.h"
#include "gin/globals.h"
#include "get_settings.h"

using namespace Rcpp;

//' Run shake.
//' 
//' @description Run shake.
//' 
//' @param X A matrix with the genotypes.
//' @param Y A vector with the phenptypes.
//' @param network A sparse matrix containing the adjacency matrix of the network.
//' @param userSettings A named list with the settings.
//' @return An object with the shake results.
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

#endif //SHAKE_TEST_ASSOCIATIONS
