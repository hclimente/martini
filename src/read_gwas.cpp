// [[Rcpp::interfaces(r,cpp)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/io/CPlinkParser.h"
#include "gin/io/CSconesIO.h"

using namespace Rcpp;

//[[Rcpp::export]]
List read_gwas(std::string pedBasename, std::string phenoFile, std::string netPath, unsigned int encoding, double maf) {
  GWASData data;
  GWASData tmpData;

  CPlinkParser::readPEDFile(pedBasename + ".ped", &tmpData);
  CPlinkParser::readMAPFile(pedBasename + ".map", &tmpData);
  CPlinkParser::readPhenotypeFile(phenoFile, &tmpData);
  CGWASDataHelper::encodeHeterozygousData(&tmpData, encoding);
  CGWASDataHelper::filterSNPsByMAF(&tmpData, maf);
  data = CGWASDataHelper::removeSamples4MissingData(tmpData, 0);

  CSconesIO::readSparseNetworkFile(netPath, &data);

  return Rcpp::List::create(Rcpp::Named("X") = data.X,
                            Rcpp::Named("Y") = data.Y.col(0),
                            Rcpp::Named("ids") = data.snp_identifiers,
                            Rcpp::Named("net") = data.network);

}
