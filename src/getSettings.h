#ifndef RSCONES2_GET_SETTINGS
#define RSCONES2_GET_SETTINGS

#include <Rcpp.h>
#include <RcppEigen.h>
#include "CEasyGWAS/gwas/CScones.h"

using namespace Rcpp;

CSconesSettings getSettings(List);

#endif //RSCONES2_GET_SETTINGS
