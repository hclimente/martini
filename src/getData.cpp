// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "CEasyGWAS/io/CPlinkParser.h"
#include "CEasyGWAS/io/CSconesIO.h"

using namespace Rcpp;

//[[Rcpp::export]]
List getData(std::string filesPath) {
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

  return Rcpp::List::create(Rcpp::Named("Y") = data.Y.col(0),
                            Rcpp::Named("X") = data.X,
                            Rcpp::Named("network") = data.network);

}
