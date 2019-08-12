#include <Rcpp.h>
using namespace Rcpp;

#ifdef __cplusplus
extern "C" {
#endif

#include "skew.h"

#ifdef __cplusplus
}
#endif

//' Call skew_eval_c function in skew.c
//'
//' @param mode A scalar, the mode of the target distribution
//' @param h A vector, the first five derivatives of the log target distribution
//' @param mu A scalar, the prior mean
//' @param omega A scalar, the prior precision
//' @param z A vector of points of evaluation
//' @return
//' A vector of evaluations of the log normalized target density
//' @export
// [[Rcpp::export]]
NumericVector skew_eval_cpp(int n_grid_points, int is_v, double mode, NumericVector h, double mu, double omega, NumericVector z)
{
  NumericVector log_f(skew_eval_c(n_grid_points, is_v, mode, h, mu, omega, z));
  return log_f;
}

//' Call skew_draw_c function in skew.c
//'
//' @param mode A scalar, the mode of the target distribution
//' @param h A vector, the first five derivatives of the log target distribution
//' @param mu A scalar, the prior mean
//' @param omega A scalar, the prior precision
//' @param n_draws A scalar integer, the number of draws to make
//' @return
//' A vector of draws from the target distribution
//' @export
// [[Rcpp::export]]
NumericVector skew_draw_cpp(int n_grid_points, int is_v, double mode, NumericVector h, double mu, double omega, int n_draws)
{
  NumericVector draws(skew_draw_c(n_grid_points, is_v, mode, h, mu, omega, n_draws));
  return draws;
}

/*** R
skew_eval_cpp(6, TRUE, 5, c(0, 0.4, 0.01, 0, 0), 4, 1.0, seq(0, 10, by=0.1))
skew_draw_cpp(6, TRUE, 5, c(0, 0.4, 0.01, 0, 0), 4, 1.0, 1000)
*/
