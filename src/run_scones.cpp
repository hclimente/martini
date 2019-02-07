#ifndef MARTINI_RUN_SCONES
#define MARTINI_RUN_SCONES

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/feature_selection/scones.h"

using namespace Rcpp;

//' Run shake.
//' 
//' @description Run scones.
//' 
//' @param c A vector with the association of each SNP with the phenotype.
//' @param eta A numeric with the value of the eta parameter.
//' @param lambda A numeric with the value of the eta parameter.
//' @param W A sparse matrix with the connectivity.
//' @return A list with vector indicating if the feature was selected and the 
//' objective score.
// [[Rcpp::export]]
Eigen::VectorXd run_scones(Eigen::VectorXd c, double eta, double lambda, 
                           Eigen::SparseMatrix<double,Eigen::ColMajor> W) {
  
  Scones s(c, eta, lambda, &W);
  s.selectSnps();
  
  return s.selected();
  
}

#endif //MARTINI_RUN_SCONES
