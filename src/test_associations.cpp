#ifndef RSCONES2_TEST_ASSOCIATIONS
#define RSCONES2_TEST_ASSOCIATIONS

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
List test_associations(int statistic, std::string filesPath, double lambda, double eta) {

  CSconesSettings settings;
  settings = CSconesSettings();

  // custom settings
  // set up specific lambda and eta
  VectorXd l(1);
  l(0) = lambda;
  VectorXd e(1);
  e(0) = eta;
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

  cout << "Testing associations.\n";
  cout << lambda << "\t" << eta << "\n";
  scones.test_associations(lambda, eta);

  cout << "Collecting results.\n";
  VectorXd indicator = scones.getIndicatorVector();
  VectorXd terms = scones.getObjectiveFunctionTerms(lambda,eta);
  double objectiveScore = scones.getObjectiveScore();
  VectorXd scores = scones.getScoreStatistic();
  SparseMatrixXd W = scones.getW();

  return Rcpp::List::create(Rcpp::Named("indicator") = indicator,
                            Rcpp::Named("terms") = terms,
                            Rcpp::Named("objective") = objectiveScore,
                            Rcpp::Named("scores") = scores,
                            Rcpp::Named("W") = W,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("eta") = eta);
}

#endif //RSCONES2_TEST_ASSOCIATIONS
