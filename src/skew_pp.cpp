#include <Rcpp.h>
using namespace Rcpp;

#ifdef __cplusplus
extern "C" {
#endif

#include "multi_skew.h"

#ifdef __cplusplus
}
#endif

//' Call skew_eval_c function in skew.c
//'
//' @param n_grid_points, number of spline grid points in f_u approximation
//' @param code, 0 for spline_skew_draw without v transformation, 1 for spline_skew_draw with v transformation,
//' 2 for old skew_draw
//' @param mode A scalar, the mode of the target distribution
//' @param h A vector, the first five derivatives of the log target distribution
//' @param mu A scalar, the prior mean
//' @param omega A scalar, the prior precision
//' @param z A vector of points of evaluation
//' @return
//' A vector of evaluations of the log normalized target density
//' @export
// [[Rcpp::export]]
NumericVector skew_eval_cpp(int n_grid_points, int code, double mode, NumericVector h, double mu, double omega, NumericVector z)
{
  NumericVector log_f(skew_eval_c(n_grid_points, code, mode, h, mu, omega, z));
  return log_f;
}

//' Call skew_draw_c function in skew.c
//'
//' @param n_grid_points, number of spline grid points in f_u approximation
//' @param code, 0 for spline_skew_draw without v transformation, 1 for spline_skew_draw with v transformation,
//' 2 for old skew_draw
//' @param mode A scalar, the mode of the target distribution
//' @param h A vector, the first five derivatives of the log target distribution
//' @param mu A scalar, the prior mean
//' @param omega A scalar, the prior precision
//' @param n_draws A scalar integer, the number of draws to make
//' @return
//' A list consisting of a vector of draws from the target distribution and a vector of corresponding
//' log density evaluations.
//' @export
// [[Rcpp::export]]
List skew_draw_cpp(int n_grid_points, int code, double mode, NumericVector h, double mu, double omega, int n_draws)
{
  List draws_and_evals(skew_draw_c(n_grid_points, code, mode, h, mu, omega, n_draws));
  return draws_and_evals;
}

