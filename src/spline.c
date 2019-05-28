#include <R.h>
#include <Rinternals.h>

#include "spline.h"

SEXP spline_eval_c(SEXP p, SEXP m, SEXP u) {
  double *p_ptr, *m_ptr, *u_ptr, *result_ptr;
  SEXP result;
  p = PROTECT(coerceVector(p, REALSXP));
  m = PROTECT(coerceVector(m, REALSXP));
  u = PROTECT(coerceVector(u, REALSXP));
  result = PROTECT(allocVector(REALSXP, length(u)));
  p_ptr = REAL(p);
  m_ptr = REAL(m);
  u_ptr = REAL(u);
  result_ptr = REAL(result);
  spline_eval(length(p), p_ptr, m_ptr, u_ptr, length(u), result_ptr);
  UNPROTECT(4);
  return result;
}

SEXP spline_draw_c(SEXP p, SEXP m, int n_draws) {
  double *p_ptr, *m_ptr, *result_ptr;
  SEXP result;
  p = PROTECT(coerceVector(p, REALSXP));
  m = PROTECT(coerceVector(m, REALSXP));
  result = PROTECT(allocVector(REALSXP, n_draws));
  p_ptr = REAL(p);
  m_ptr = REAL(m);
  result_ptr = REAL(result);
  spline_draw(length(p), p_ptr, m_ptr, n_draws, result_ptr);
  UNPROTECT(3);
  return result;
}

