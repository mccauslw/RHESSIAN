#ifndef MULTI_SKEW
#define MULTI_SKEW

// For use in C++ wrapper used for testing
SEXP skew_eval_c(int grid_points, int code, double mode, SEXP h, double mu, double omega, SEXP z);
SEXP skew_draw_c(int grid_points, int code, double mode, SEXP h, double mu, double omega, int n_draws);

#endif /* MULTI_SKEW */
