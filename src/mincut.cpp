#ifndef MARTINI_MINCUTC
#define MARTINI_MINCUTC

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "maxflow/maxflow.h"

using namespace Rcpp;

//' Maxflow algorithm
//' 
//' @description Run the maxflow algorithm.
//' @param A A sparse matrix with the connectivity.
//' @param As A vector containing the edges to the source.
//' @param At A vector containing the edges to the sink.
//' @return A list with vector indicating if the feature was selected and the 
//' objective score.
// [[Rcpp::export]]
LogicalVector maxflow(Eigen::SparseMatrix<double,Eigen::ColMajor> const &A,
                      Eigen::VectorXd const &As,
                      Eigen::VectorXd const &At) {
  
  LogicalVector selected(A.rows());
  
  // create graph out of adjacency matrix
  typedef Graph<double, double, double> MaxGraph;
  MaxGraph *g = new MaxGraph(A.rows(), A.nonZeros());
  
  // initialize nodes
  g->add_node(selected.length());
  
  // add edge weights from the original graph A
  for(int k = 0; k < A.outerSize(); k++)
    for(Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
      g->add_edge(it.row(), it.col(), it.value(), 0.0);
  
  // add edges to the source
  for(int k = 0; k < A.rows(); k++)
    if(As(k) != 0)
      g->add_tweights(k, As(k), 0.0);
    
  // add edges to the sink
  for(int k = 0; k < A.rows(); k++)
    if(At(k) != 0)
      g->add_tweights(k, 0.0, At(k));
  
  // run maxflow algorithm
  g->maxflow();
  
  // create indicator_vector
  for(int i = 0; i < selected.length(); i++)
    selected[i] = g->what_segment(i);
  
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
LogicalVector mincut_c(Eigen::VectorXd c, double eta, double lambda, 
                       Eigen::SparseMatrix<double,Eigen::ColMajor> W) {
  
  Eigen::SparseMatrix<double,Eigen::ColMajor> A = lambda * W;
  Eigen::VectorXd c_t = c.array() - eta;
  
  // add source and sink
  // connect negative c values to source
  Eigen::VectorXd As = (c_t.array() > 0).select(0, -c_t);
  // connect positive c values to sink
  Eigen::VectorXd At = (c_t.array() <= 0).select(0, c_t);
  
  // compute maxflow
  LogicalVector selected = maxflow(A, As, At);
  return(selected);
  
}

#endif //MARTINI_MINCUTC
