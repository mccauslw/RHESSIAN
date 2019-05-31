#ifndef SKEWGRID
#define SKEWGRID

// Structure to hold grid values that never change
typedef struct {
  int K;
  double *u;
  double *v;
  double *z;
  double *f_v;
  double *f_v_prime;
  double *z_plus;
  double *z_minus;
} Grid;

#endif // SKEWGRID
