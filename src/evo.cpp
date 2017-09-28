#ifndef MARTINI_EVO
#define MARTINI_EVO

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/gwas/CGWASData.h"
#include "gin/feature_selection/shake.h"
#include "gin/globals.h"
#include "gin/settings.h"
using namespace Rcpp;

//' Run evo.
//' 
//' @description Run evo.
//' 
//' @param X A matrix with the genotypes.
//' @param Y A vector with the phenptypes.
//' @param network A sparse matrix containing the adjacency matrix of the network.
//' @param userSettings A named list with the settings.
//' @return An object with the evo results.
// [[Rcpp::export]]
Rcpp::List evo(Eigen::MatrixXd X, Eigen::VectorXd Y, Eigen::SparseMatrix<double,Eigen::ColMajor> network, Rcpp::List userSettings) {

  VectorXd etas(Rcpp::as<Eigen::VectorXd>(userSettings["etas"]));
  VectorXd lambdas(Rcpp::as<Eigen::VectorXd>(userSettings["lambdas"]));
  
  Settings s("", "", userSettings["encoding"], userSettings["modelScore"], userSettings["associationScore"], etas, lambdas, "");
  
  Shake experiment(X, Y, network);
  experiment.setDebug(userSettings["debug"]);
  
  if (s.etas().rows() == 0) {
    experiment.selectHyperparameters(10, s.modelScore(), s.associationScore());
  } else {
    experiment.selectHyperparameters(10, s.modelScore(), s.associationScore(), s.etas(), s.lambdas());
  }
  
  experiment.selectSNPs();
  
  return Rcpp::List::create(Rcpp::Named("selected") = experiment.selectedSnps(),
                            Rcpp::Named("c") = experiment.c(),
                            Rcpp::Named("lambda") = experiment.bestLambda(),
                            Rcpp::Named("eta") = experiment.bestEta());
  
}

#endif //MARTINI_EVO
