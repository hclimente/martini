// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "CEasyGWAS/gwas/CScones.h"

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd maxflow(Eigen::SparseMatrix<double,Eigen::ColMajor> lW, Eigen::MatrixXd A) {

  CScones scones;

  VectorXd indicator_vector;

  // Matrix A containing As and At
  MatrixXd A(__n_features,2);
  //connect positive c values to sink
  VectorXd pos_c = (c.array()<=0).select(0,c);
  //connect negative c values to source
  VectorXd neg_c = -c;
  neg_c = (neg_c.array()<=0).select(0,neg_c);
  //Store data
  A.col(0) = neg_c;
  A.col(1) = pos_c;

  //compute maxflow
  scones.maxflow(lW, A, &indicator_vector);

  return indicator_vector;
}
