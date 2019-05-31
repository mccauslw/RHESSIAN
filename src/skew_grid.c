#include "skew_grid.h"

// These depend only on n_knots and use_delta
static double u[6];
static double v[6];
static double z[6];
static double f_v[6];
static double f_v_prime[6];
static double z_plus[36];
static double z_minus[36];

Grid g6 = {5, u, v, z, f_v, f_v_prime, z_plus, z_minus};
