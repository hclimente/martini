#ifndef MARTINI_MINCUTC
#define MARTINI_MINCUTC

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "maxflow/maxflow.h"

using namespace Rcpp;

// [[Rcpp::export]]
Eigen::VectorXd maxflow(Eigen::MatrixXd const &A,
                        Eigen::SparseMatrix<double,Eigen::ColMajor> const &W) {
  
  Eigen::VectorXd selected = Eigen::VectorXd::Ones(W.rows());
  
  // create graph out of adjacency matrix
  typedef Graph<double, double, double> MaxGraph;
  MaxGraph *g = new MaxGraph(W.rows(), W.nonZeros());
  
  // initialize nodes
  g->add_node(selected.rows());
  
  //traverse the sparse adjacency matrix A
  for(long long k=0; k<W.outerSize(); k++) {
    for(Eigen::SparseMatrix<double>::InnerIterator it(W,k); it; ++it) {
      g->add_edge(it.row(), it.col(), it.value(), 0.0);
    }
  }
  
  // traverse the T matrix
  for(long long i = 0; i < 2; i++) {
    for(long long k = 0; k < A.rows(); k++) {
      if(i==0) {
        if(A(k,i) != 0) {
          g->add_tweights(k, A(k,i), 0.0);
        }
      } else {
        if(A(k,i) != 0) {
          g->add_tweights(k, 0.0, A(k,i));
        }
      }
    }
  }
  
  // run maxflow algorithm
  g->maxflow();
  
  // create indicator_vector
  for(long long i = 0; i < selected.rows(); i++)
    selected(i) = g->what_segment(i);
  
  // delete graph
  delete g;
  
  return(selected);
}

//' Min-cut algorithm
//' 
//' @description Run the mincut algorithm.
//' @param c A vector with the association of each SNP with the phenotype.
//' @param eta A numeric with the value of the eta parameter.
//' @param lambda A numeric with the value of the eta parameter.
//' @param W A sparse matrix with the connectivity.
//' @return A list with vector indicating if the feature was selected and the 
//' objective score.
// [[Rcpp::export]]
Eigen::VectorXd mincut_c(Eigen::VectorXd c, double eta, double lambda, 
                         Eigen::SparseMatrix<double,Eigen::ColMajor> W) {
  
  W = lambda * W;
  long n_features = c.rows();
  
  Eigen::VectorXd c_t = c.array() - eta;
  
  // Add source and sink
  // Matrix A containing As and At
  Eigen::MatrixXd A(n_features, 2);
  // connect positive c values to sink
  Eigen::VectorXd pos_c = (c_t.array() <= 0).select(0, c_t);
  // connect negative c values to source
  Eigen::VectorXd neg_c = - c_t;
  neg_c = (c_t.array() > 0).select(0, neg_c);
  //Store data
  A.col(0) = neg_c;
  A.col(1) = pos_c;
  
  //compute maxflow
  Eigen::VectorXd selected = maxflow(A, W);
  return(selected);
  
}

#endif //MARTINI_MINCUTC
