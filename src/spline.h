#ifndef SPLINE
#define SPLINE

// For internal C use

void spline_eval(int n_knots, double *p, double *m, int n_evals, double *u, double *f_u);
void spline_draw(int n_knots, double *p, double *m, int n_draws, double *u);

// For use in C++ wrapper spline_pp.cpp used for testing
SEXP spline_eval_c(SEXP p, SEXP m, SEXP u);
SEXP spline_draw_c(SEXP p, SEXP m, int n_draws);

#endif /* SPLINE */
