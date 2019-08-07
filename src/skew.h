#ifndef SKEW
#define SKEW

#include "skew_grid.h"

// For internal C use
void skew_draw_eval(int is_draw, int n_draws, Grid *g,
                    double mode, double *h, double mu, double omega,
                    double *z, double *ln_f);

// For use in C++ wrapper used for testing
SEXP skew_eval_c(double mode, SEXP h, double mu, double omega, SEXP z);
SEXP skew_draw_c(double mode, SEXP h, double mu, double omega, int n_draws);

#endif /* SKEW */
