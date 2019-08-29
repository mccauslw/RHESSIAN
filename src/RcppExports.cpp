// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// alias_cpp
NumericVector alias_cpp(NumericVector x, int n_draws);
RcppExport SEXP _RHESSIAN_alias_cpp(SEXP xSEXP, SEXP n_drawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n_draws(n_drawsSEXP);
    rcpp_result_gen = Rcpp::wrap(alias_cpp(x, n_draws));
    return rcpp_result_gen;
END_RCPP
}
// skew_eval_cpp
NumericVector skew_eval_cpp(int n_grid_points, int code, double mode, NumericVector h, double mu, double omega, NumericVector z);
RcppExport SEXP _RHESSIAN_skew_eval_cpp(SEXP n_grid_pointsSEXP, SEXP codeSEXP, SEXP modeSEXP, SEXP hSEXP, SEXP muSEXP, SEXP omegaSEXP, SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_grid_points(n_grid_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type code(codeSEXP);
    Rcpp::traits::input_parameter< double >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(skew_eval_cpp(n_grid_points, code, mode, h, mu, omega, z));
    return rcpp_result_gen;
END_RCPP
}
// skew_draw_cpp
List skew_draw_cpp(int n_grid_points, int code, double mode, NumericVector h, double mu, double omega, int n_draws);
RcppExport SEXP _RHESSIAN_skew_draw_cpp(SEXP n_grid_pointsSEXP, SEXP codeSEXP, SEXP modeSEXP, SEXP hSEXP, SEXP muSEXP, SEXP omegaSEXP, SEXP n_drawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_grid_points(n_grid_pointsSEXP);
    Rcpp::traits::input_parameter< int >::type code(codeSEXP);
    Rcpp::traits::input_parameter< double >::type mode(modeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< int >::type n_draws(n_drawsSEXP);
    rcpp_result_gen = Rcpp::wrap(skew_draw_cpp(n_grid_points, code, mode, h, mu, omega, n_draws));
    return rcpp_result_gen;
END_RCPP
}
// spline_eval_cpp
NumericVector spline_eval_cpp(NumericVector p, NumericVector m, NumericVector u);
RcppExport SEXP _RHESSIAN_spline_eval_cpp(SEXP pSEXP, SEXP mSEXP, SEXP uSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    rcpp_result_gen = Rcpp::wrap(spline_eval_cpp(p, m, u));
    return rcpp_result_gen;
END_RCPP
}
// spline_draw_cpp
NumericVector spline_draw_cpp(NumericVector p, NumericVector m, int n_draws);
RcppExport SEXP _RHESSIAN_spline_draw_cpp(SEXP pSEXP, SEXP mSEXP, SEXP n_drawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n_draws(n_drawsSEXP);
    rcpp_result_gen = Rcpp::wrap(spline_draw_cpp(p, m, n_draws));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RHESSIAN_alias_cpp", (DL_FUNC) &_RHESSIAN_alias_cpp, 2},
    {"_RHESSIAN_skew_eval_cpp", (DL_FUNC) &_RHESSIAN_skew_eval_cpp, 7},
    {"_RHESSIAN_skew_draw_cpp", (DL_FUNC) &_RHESSIAN_skew_draw_cpp, 7},
    {"_RHESSIAN_spline_eval_cpp", (DL_FUNC) &_RHESSIAN_spline_eval_cpp, 3},
    {"_RHESSIAN_spline_draw_cpp", (DL_FUNC) &_RHESSIAN_spline_draw_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RHESSIAN(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
