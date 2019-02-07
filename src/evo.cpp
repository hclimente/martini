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
//' @param W A sparse matrix containing the adjacency matrix of the network.
//' @param opts A named list with the settings.
//' @return An object with the evo results.
// [[Rcpp::export]]
Rcpp::List evo(Eigen::MatrixXd X, Eigen::VectorXd Y, Eigen::SparseMatrix<double,
               Eigen::ColMajor> W, Rcpp::List opts) {

  VectorXd etas(Rcpp::as<Eigen::VectorXd>(opts["etas"]));
  VectorXd lambdas(Rcpp::as<Eigen::VectorXd>(opts["lambdas"]));
  
  Settings s("", "", 0, opts["modelScore"], 
             opts["associationScore"], etas, lambdas, "");
  
  Shake exp(&X, &Y, &W);
  exp.setDebug(opts["debug"]);
  
  if (s.etas().rows() == 0) {
    exp.selectHyperparameters(10, s.modelScore(), s.associationScore());
  } else {
    exp.selectHyperparameters(10, s.modelScore(), s.associationScore(), 
                              s.etas(), s.lambdas());
  }
  
  exp.selectSNPs();
  
  return Rcpp::List::create(Rcpp::Named("selected") = exp.selectedSnps(),
                            Rcpp::Named("c") = exp.c(),
                            Rcpp::Named("lambdas") = exp.grid()->lambdas(),
                            Rcpp::Named("etas") = exp.grid()->etas(),
                            Rcpp::Named("lambda") = exp.bestLambda(),
                            Rcpp::Named("eta") = exp.bestEta(),
                            Rcpp::Named("grid") = exp.grid()->scoredFolds());
  
}

#endif //MARTINI_EVO
