#ifndef RSCONES2_TEST_ASSOCIATIONS_GRIDSEARCH
#define RSCONES2_TEST_ASSOCIATIONS_GRIDSEARCH

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
List testAssociationsGridsearch(int statistic, std::string filesPath, unsigned int gridparams, int griddepth, unsigned int criterion){

  CSconesSettings settings;
  settings = CSconesSettings();

  // custom settings
  settings.test_statistic = statistic;
  settings.selection_criterion = criterion;
  settings.nParameters = gridparams;
  settings.gridsearch_depth = griddepth;
  settings.lambdas = VectorXd::Zero(settings.nParameters);
  settings.etas = VectorXd::Zero(settings.nParameters);

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
  scones.test_associations();

  double eta_f = scones.getBestEta();
  double lambda_f =  scones.getBestLambda();
  VectorXd indicator = scones.getIndicatorVector();
  VectorXd terms = scones.getObjectiveFunctionTerms(lambda_f, eta_f);
  // double objectiveScore = scones.getObjectiveScore();
  VectorXd scores = scones.getScoreStatistic();
  SparseMatrixXd W = scones.getW();

  return Rcpp::List::create(Rcpp::Named("indicator") = indicator,
                            Rcpp::Named("terms") = terms,
                            // Rcpp::Named("objective") = objectiveScore,
                            Rcpp::Named("scores") = scores,
                            Rcpp::Named("W") = W,
                            Rcpp::Named("eta") = eta_f,
                            Rcpp::Named("lambda") = lambda_f);

}

#endif //RSCONES2_TEST_ASSOCIATIONS_GRIDSEARCH
