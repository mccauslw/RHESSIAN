#ifndef SKEW
#define SKEW

#include "skew_grid.h"

// For internal C use
void skew_draw_eval(double mode, double *h, double mu, double omega,
                    int n, double *z, double *ln_f,
                    int is_draw, Grid *g);

// For use in C++ wrapper used for testing
SEXP skew_eval_c(double mode, SEXP h, double mu, double omega, SEXP z);
SEXP skew_draw_c(double mode, SEXP h, double mu, double omega, int n_draws);

#endif /* SKEW */
