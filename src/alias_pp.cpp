#include <Rcpp.h>
using namespace Rcpp;

#ifdef __cplusplus
extern "C" {
#endif

#include "alias.h"

#ifdef __cplusplus
}
#endif

//' Call alias_c function from C
//'
//' @param x A vector of unnormalized probabilities
//' @param n_draws Integer number of draws from the discrete distribution
//' @return
//' An integer-valued vector of draws from the discrete distribution
//' @export
// [[Rcpp::export]]
NumericVector alias_cpp(NumericVector x, int n_draws) {
  NumericVector result(alias_c(x, n_draws));
  return result;
}

/*** R
alias_cpp(c(1,2,3,4,5,4,3,2,1), 100)
*/
