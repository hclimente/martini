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
List test_associations(int statistic, std::string filesPath, double eta, double lambda) {

  CSconesSettings settings;
  CScones scones;
  GWASData tmpData;

  GWASData data;

  string genotype_str = filesPath + "genotype";
  string phenotype_str = filesPath + "phenotype.txt";
  string network_str = filesPath + "network.txt";
  uint encoding = 0;
  float64 maf = 0.05;

  CPlinkParser::readPEDFile(genotype_str + ".ped", &data);
  CPlinkParser::readMAPFile(genotype_str + ".map", &data);
  CPlinkParser::readPhenotypeFile(phenotype_str,&data);
  CGWASDataHelper::encodeHeterozygousData(&data,encoding);
  CGWASDataHelper::filterSNPsByMAF(&data,maf);
  CSconesIO::readSparseNetworkFile(network_str,&data);
  tmpData = CGWASDataHelper::removeSamples4MissingData(data,0);

  settings = CSconesSettings();
  settings.folds = 10;
  settings.seed = 0;
  settings.selection_criterion = CONSISTENCY;
  settings.selection_ratio = 0.8;
  settings.test_statistic = statistic;
  settings.nParameters = 10;
  settings.evaluateObjective = false;
  settings.dump_intermediate_results = true;
  settings.dump_path = "tmp/";

  // set up specific lambda and eta
  VectorXd l(1);
  l(0) = lambda;
  VectorXd e(1);
  e(0) = eta;
  settings.lambdas = l;
  settings.etas = e;
  // avoid gridsearch
  settings.autoParameters = false;

  scones = CScones(tmpData.Y.col(0),tmpData.X,tmpData.network, settings);

  scones.test_associations(eta, lambda);

  VectorXd indicator = scones.getIndicatorVector();
  VectorXd terms = scones.getObjectiveFunctionTerms(lambda,eta);
  VectorXd scores = scones.getScoreStatistic();
  SparseMatrixXd W = scones.getW();

  return Rcpp::List::create(Rcpp::Named("indicator") = indicator,
                            Rcpp::Named("terms") = terms,
                            Rcpp::Named("scores") = scores,
                            Rcpp::Named("W") = W,
                            Rcpp::Named("eta") = eta,
                            Rcpp::Named("lambda") = lambda,
                            Rcpp::Named("eta") = eta);
}
