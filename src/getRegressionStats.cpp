#ifndef RSCONES2_GET_REGRESSION_STATS
#define RSCONES2_GET_REGRESSION_STATS

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "CEasyGWAS/regression/CRegression.h"
#include "CEasyGWAS/globals.h"
#include "getSettings.h"

using namespace Rcpp;

// [[Rcpp::export]]
List getRegressionStats(Eigen::MatrixXd X, Eigen::VectorXd Y) {

  double BIC;
  double AIC;
  double AICc;
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
    BIC = logistic_regression.getBIC();
    AIC = logistic_regression.getAIC();
    AICc = logistic_regression.getAICc();
    logLikelihood = logistic_regression.getLogLikelihood();
  } else { //if phenotype is continuous select LinearRegression
    CLinearRegression linear_regression;
    linear_regression = CLinearRegression();
    linear_regression.fit(Y, X);
    BIC = linear_regression.getBIC();
    AIC = linear_regression.getAIC();
    AICc = linear_regression.getAICc();
    logLikelihood = linear_regression.getLogLikelihood();
  }

  return Rcpp::List::create(Rcpp::Named("BIC") = BIC,
                            Rcpp::Named("AIC") = AIC,
                            Rcpp::Named("AICc") = AICc,
                            Rcpp::Named("logLikelihood") = logLikelihood);

}

#endif //RSCONES2_GET_REGRESSION_STATS
