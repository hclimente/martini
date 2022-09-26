// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/martini.h"
#include <RcppEigen.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// maxflow
LogicalVector maxflow(Eigen::SparseMatrix<double,Eigen::ColMajor> const& A, Eigen::VectorXd const& As, Eigen::VectorXd const& At);
static SEXP _martini_maxflow_try(SEXP ASEXP, SEXP AsSEXP, SEXP AtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double,Eigen::ColMajor> const& >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd const& >::type As(AsSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd const& >::type At(AtSEXP);
    rcpp_result_gen = Rcpp::wrap(maxflow(A, As, At));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _martini_maxflow(SEXP ASEXP, SEXP AsSEXP, SEXP AtSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_martini_maxflow_try(ASEXP, AsSEXP, AtSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}
// mincut_c
LogicalVector mincut_c(Eigen::VectorXd c, double eta, double lambda, Eigen::SparseMatrix<double,Eigen::ColMajor> W);
static SEXP _martini_mincut_c_try(SEXP cSEXP, SEXP etaSEXP, SEXP lambdaSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type c(cSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double,Eigen::ColMajor> >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(mincut_c(c, eta, lambda, W));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _martini_mincut_c(SEXP cSEXP, SEXP etaSEXP, SEXP lambdaSEXP, SEXP WSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_martini_mincut_c_try(cSEXP, etaSEXP, lambdaSEXP, WSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    bool rcpp_isLongjump_gen = Rcpp::internal::isLongjumpSentinel(rcpp_result_gen);
    if (rcpp_isLongjump_gen) {
        Rcpp::internal::resumeJump(rcpp_result_gen);
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _martini_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("LogicalVector(*maxflow)(Eigen::SparseMatrix<double,Eigen::ColMajor> const&,Eigen::VectorXd const&,Eigen::VectorXd const&)");
        signatures.insert("LogicalVector(*mincut_c)(Eigen::VectorXd,double,double,Eigen::SparseMatrix<double,Eigen::ColMajor>)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _martini_RcppExport_registerCCallable() { 
    R_RegisterCCallable("martini", "_martini_maxflow", (DL_FUNC)_martini_maxflow_try);
    R_RegisterCCallable("martini", "_martini_mincut_c", (DL_FUNC)_martini_mincut_c_try);
    R_RegisterCCallable("martini", "_martini_RcppExport_validate", (DL_FUNC)_martini_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_martini_maxflow", (DL_FUNC) &_martini_maxflow, 3},
    {"_martini_mincut_c", (DL_FUNC) &_martini_mincut_c, 4},
    {"_martini_RcppExport_registerCCallable", (DL_FUNC) &_martini_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_martini(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
