#ifndef SKEW
#define SKEW

// For internal C use
void skew_draw_eval(double mode, double *h, double mu, double omega,
                    int n, double *z, double *ln_f, int is_draw, int n_knots);

// For use in C++ wrapper used for testing
SEXP skew_eval_c(double mode, SEXP h, double mu, double omega, SEXP z);
SEXP skew_draw_c(double mode, SEXP h, double mu, double omega, int n_draws);

#endif /* SKEW */
