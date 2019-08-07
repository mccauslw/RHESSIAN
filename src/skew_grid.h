#ifndef SKEWGRID
#define SKEWGRID

// Structure to hold grid values that never change
typedef struct {
  int K;
  double *u;
  double *v;
  double *x;
  double *c;
  double *f_v;
  double *f_v_prime;
  double *x_plus;
  double *x_minus;
} Grid;

#endif // SKEWGRID
