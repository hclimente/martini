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
  settings.nParameters = 0;
  settings.evaluateObjective = true;
  settings.dump_intermediate_results = true;
  settings.dump_path = "tmp/";

  // set up specific lambda and eta
  VectorXd l(1);
  l(0) = 0;
  VectorXd e(1);
  e(0) = 0;
  settings.lambdas = l;
  settings.etas = e;
  // avoid gridsearch
  settings.autoParameters = false;

  scones = CScones(tmpData.Y.col(0),tmpData.X,tmpData.network, settings);

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
