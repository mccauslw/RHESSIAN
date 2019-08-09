#ifndef SKEW
#define SKEW

#include "skew_grid.h"

// For internal C use
void skew_draw_eval(int is_draw, int n_draws, int n_grid_points, int is_v,
                    double mode, double *h, double mu, double omega,
                    double *z, double *ln_f);

// For use in C++ wrapper used for testing
SEXP skew_eval_c(int grid_points, double mode, SEXP h, double mu, double omega, SEXP z);
SEXP skew_draw_c(int grid_points, double mode, SEXP h, double mu, double omega, int n_draws);

#endif /* SKEW */
