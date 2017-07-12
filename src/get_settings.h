#ifndef RSCONES2_GET_SETTINGS
#define RSCONES2_GET_SETTINGS

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/gwas/CScones.h"

using namespace Rcpp;

CSconesSettings get_settings(List);

#endif //RSCONES2_GET_SETTINGS
