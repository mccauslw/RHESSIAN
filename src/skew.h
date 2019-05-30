#ifndef SKEW
#define SKEW

// Structure to hold grid values that never change
typedef struct {
  int n_knots;
  double *u;
  double *v;
  double *z;
  double *f_v;
  double *f_v_prime;
  double *z_plus;
  double *z_minus;
} Grid;

// For internal C use
void skew_draw_eval(double mode, double *h, double mu, double omega,
                    int n, double *z, double *ln_f, int is_draw, int n_knots);

// For use in C++ wrapper used for testing
SEXP skew_eval_c(double mode, SEXP h, double mu, double omega, SEXP z);
SEXP skew_draw_c(double mode, SEXP h, double mu, double omega, int n_draws);

#endif /* SKEW */
