#ifndef RSCONES2_GET_REGRESSION_STATS
#define RSCONES2_GET_REGRESSION_STATS

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/regression/CRegression.h"
#include "gin/globals.h"
#include "get_settings.h"

using namespace Rcpp;

// [[Rcpp::export]]
List get_regression_stats(Eigen::MatrixXd X, Eigen::VectorXd Y) {

  double dBIC;
  double dAIC;
  double dAICc;
  double logLikelihood;

  bool binary = true;
  for(int64 i=0; i<Y.rows();i++) {
    if(!(Y(i)==0 || Y(i)==1)) {
      binary = false;
      break;
    }
  }

  if(binary) { //if phenotype is binary select LogisticRegression
    CLogisticRegression logistic_regression;
    logistic_regression = CLogisticRegression();
    logistic_regression.fit(Y, X);
    dBIC = logistic_regression.getBIC();
    dAIC = logistic_regression.getAIC();
    dAICc = logistic_regression.getAICc();
    logLikelihood = logistic_regression.getLogLikelihood();
  } else { //if phenotype is continuous select LinearRegression
    CLinearRegression linear_regression;
    linear_regression = CLinearRegression();
    linear_regression.fit(Y, X);
    dBIC = linear_regression.getBIC();
    dAIC = linear_regression.getAIC();
    dAICc = linear_regression.getAICc();
    logLikelihood = linear_regression.getLogLikelihood();
  }

  return Rcpp::List::create(Rcpp::Named("BIC") = dBIC,
                            Rcpp::Named("AIC") = dAIC,
                            Rcpp::Named("AICc") = dAICc,
                            Rcpp::Named("logLikelihood") = logLikelihood);

}

#endif //RSCONES2_GET_REGRESSION_STATS
