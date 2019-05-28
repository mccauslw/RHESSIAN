#include <Rcpp.h>
using namespace Rcpp;

#ifdef __cplusplus
extern "C" {
#endif

#include "spline.h"

#ifdef __cplusplus
}
#endif

//' Call C language spline_eval_c function
//'
//' @param p A vector of levels on grid
//' @param m A vector of derivatives on grid
//' @param u A vector of values at which to evaluate a spline density
//' @return
//' The vector of spline density values
//' @export
// [[Rcpp::export]]
NumericVector spline_eval_cpp(NumericVector p, NumericVector m, NumericVector u) {
  NumericVector result(spline_eval_c(p, m, u));
  return result;
}

//' Call C language spline_draw_c function
//'
//' @param p A vector of levels on grid
//' @param m A vector of derivatives on grid
//' @param n_draws The number of draws to make from the spline density
//' @return
//' The vector of spline density draws
//' @export
//' [[Rcpp::export]]
NumericVector spline_draw_cpp(NumericVector p, NumericVector m, int n_draws) {
  NumericVector result(spline_draw_c(p, m, n_draws));
  return result;
}

/*** R
spline_eval_cpp(c(1, 2, 1), c(1, 0, -1), seq(0, 1, by=0.01))
spline_draw_cpp(c(1, 2, 1), c(1, 0, -1), 100)
*/
