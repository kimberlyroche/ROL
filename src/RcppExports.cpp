// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// Riemann_dist_samples
Eigen::MatrixXd Riemann_dist_samples(Eigen::MatrixXd samples, int n_indiv, int n_samples_per);
RcppExport SEXP _ROL_Riemann_dist_samples(SEXP samplesSEXP, SEXP n_indivSEXP, SEXP n_samples_perSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type samples(samplesSEXP);
    Rcpp::traits::input_parameter< int >::type n_indiv(n_indivSEXP);
    Rcpp::traits::input_parameter< int >::type n_samples_per(n_samples_perSEXP);
    rcpp_result_gen = Rcpp::wrap(Riemann_dist_samples(samples, n_indiv, n_samples_per));
    return rcpp_result_gen;
END_RCPP
}
// Riemann_dist_pair
double Riemann_dist_pair(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _ROL_Riemann_dist_pair(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(Riemann_dist_pair(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ROL_Riemann_dist_samples", (DL_FUNC) &_ROL_Riemann_dist_samples, 3},
    {"_ROL_Riemann_dist_pair", (DL_FUNC) &_ROL_Riemann_dist_pair, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ROL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
