#ifndef RSCONES2_TEST_ASSOCIATIONS_CUSTOM_GRIDSEARCH
#define RSCONES2_TEST_ASSOCIATIONS_CUSTOM_GRIDSEARCH

// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "CEasyGWAS/gwas/CScones.h"
#include "CEasyGWAS/io/CSconesIO.h"
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/globals.h"

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<Rcpp::List> test_associations_custom_gridsearch(int statistic, std::string filesPath, int min, int max){

  CSconesSettings settings;
  settings = CSconesSettings();

  // custom settings
  settings.test_statistic = statistic;
  // set up specific lambda and eta
  VectorXd l(1);
  l(0) = 0;
  VectorXd e(1);
  e(0) = 0;
  settings.lambdas = l;
  settings.etas = e;
  // avoid gridsearch
  settings.autoParameters = false;

  GWASData data;
  GWASData tmpData;

  string genotype_str = filesPath + "genotype";
  string phenotype_str = filesPath + "phenotype.txt";
  string network_str = filesPath + "network.txt";
  uint encoding = 0;
  float64 maf = 0.05;

  CPlinkParser::readPEDFile(genotype_str + ".ped", &tmpData);
  CPlinkParser::readMAPFile(genotype_str + ".map", &tmpData);
  CPlinkParser::readPhenotypeFile(phenotype_str, &tmpData);
  CGWASDataHelper::encodeHeterozygousData(&tmpData, encoding);
  CGWASDataHelper::filterSNPsByMAF(&tmpData, maf);
  CSconesIO::readSparseNetworkFile(network_str, &tmpData);
  data = CGWASDataHelper::removeSamples4MissingData(tmpData, 0);

  CScones scones;
  scones = CScones(data.Y.col(0), data.X, data.network, settings);

  vector<Rcpp::List> grid;
  for(int e = min; e <= max; e++){
    for(int l = min; l <= max; l++){

      double lambda = pow(10, l);
      double eta = pow(10, e);

      scones.test_associations(lambda, eta);

      VectorXd indicator = scones.getIndicatorVector();
      VectorXd terms = scones.getObjectiveFunctionTerms(lambda, eta);
      double objectiveScore = scones.getObjectiveScore();
      VectorXd scores = scones.getScoreStatistic();
      SparseMatrixXd W = scones.getW();

      grid.push_back(Rcpp::List::create(Rcpp::Named("indicator") = indicator,
                                        Rcpp::Named("terms") = terms,
                                        Rcpp::Named("objective") = objectiveScore,
                                        Rcpp::Named("scores") = scores,
                                        Rcpp::Named("W") = W,
                                        Rcpp::Named("lambda") = lambda,
                                        Rcpp::Named("eta") = eta));

    }
  }

  return grid;
}

#endif //RSCONES2_TEST_ASSOCIATIONS_CUSTOM_GRIDSEARCH
