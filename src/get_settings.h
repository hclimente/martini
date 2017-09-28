#ifndef MARTINI_SETTINGS
#define MARTINI_SETTINGS

#include <Rcpp.h>
#include <RcppEigen.h>
#include "gin/gwas/CScones.h"

using namespace Rcpp;

CSconesSettings get_settings(List);

#endif //MARTINI_SETTINGS
