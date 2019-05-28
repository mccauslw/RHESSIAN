#include <stdlib.h>
#include <math.h>
#include "new_skew.h"
#include "alias.h"
#include "RNG.h"
#include "Phi.h"

// Used options(digits=22), log(sqrt(2/pi)) in R
const double log_root_2_by_pi = -0.2257913526447273833142;
const double delta = 0.5;
const double delta_dual = 2.0;

// Code for new skew draw.
// Documentation is in "Flexible even densities", file name flexeven.Rmd

/*
  The Grid structure stores many functions evaluated on a grid of points
  u_0, u_1, ..., u_n. See new_skew.h for details
*/
static Grid grid5_delta;
static Grid grid8;

#define n_powers 5
static double sigma[n_powers+1], a[n_powers+1];

static double f_u_0 = delta_dual;
static double fp_u_0 = -delta * delta_dual * delta_dual;

// Functions related to F_v(v) = 1-(1-v)^delta, for delta=0.5
static inline double inverse_F_v(double u) {return 1.0-(1.0-u)*(1.0-u);}
static inline double F_v(double v) {return 1.0-sqrt(1.0-v);}
static inline double f_v(double v) {return 0.5/sqrt(1.0-v);}
static inline double ln_f_v(double v) {return log(0.5) - 0.5*log(1.0-v);}

// This function sets values only depending on n, the number of grid
// points, and the boolean value use_delta, which determines whether or
// not the F_v transformation is actif or not.
static void grid_initialize(int n, int use_delta, Grid *grid)
{
	int knot, n_knots = n+1, power;
	grid->n_bins = n;
	grid->n_knots = n_knots;
	grid->use_delta = use_delta;
	grid->u = (double *) malloc(n_knots);
	grid->v = (double *) malloc(n_knots);
	grid->z = (double *) malloc(n_knots);

	grid->f_v = (double *) malloc(n_knots);
	grid->f_v_prime = (double *) malloc(n_knots);

	grid->p = (double *) malloc(n_knots);
	grid->m = (double *) malloc(n_knots);

	grid->phi_plus = (double *) malloc(n_knots);
	grid->phi_minus = (double *) malloc(n_knots);
	grid->phi_plus_prime = (double *) malloc(n_knots);
	grid->phi_minus_prime = (double *) malloc(n_knots);

	grid->z_plus = (double *) malloc(n_knots * (n_powers+1));
	grid->z_minus = (double *) malloc(n_knots * (n_powers+1));

	for (knot=0; knot<n_knots; knot++) { // WJM: case z=infty
		// Precompute values of u, v and z at current knot.
		if (knot==n_knots-1)
			grid->u[knot] = 1.0 - 0.5 / n;
		else
			grid->u[knot] = knot * 1.0 / n;
		grid->v[knot] = inverse_F_v(grid->u[knot]);
		grid->z[knot] = inverse_Phi(0.5 + 0.5*grid->v[knot]);
		grid->f_v[knot] = f_v(grid->v[knot]);
		grid->f_v_prime[knot] = 0.5*grid->f_v[knot]/(1-grid->v[knot]);

		// Precompute powers of z and -z at current knot.
		grid->z_plus[knot] = grid->z_minus[knot] = 1.0;
		for (power=1; power<=n_powers; power++) {
			grid->z_plus[knot + power*n_knots]
				= grid->z_plus[knot + (power-1)*n_knots] * grid->z[knot] / power;
			grid->z_minus[knot + power*n_knots]
				= grid->z_minus[knot + (power-1)*n_knots] * -grid->z[knot] / power;
		}
	}
}

static void grid_free(Grid *grid)
{
	free(grid->u);
	free(grid->v);
	free(grid->z);
	free(grid->f_v);
	free(grid->f_v_prime);
	free(grid->p);
	free(grid->m);
	free(grid->phi_plus);
	free(grid->phi_minus);
	free(grid->phi_plus_prime);
	free(grid->phi_minus_prime);
	free(grid->z_plus);
	free(grid->z_minus);
}

static double spline_eval(double u, Grid *grid) {
	int n = grid->n_knots;
	int i = floor(u*n);	     // Subinterval index in {0,1,...,n-1}
	double t = u*n - i;      // index of u in subinterval [i/n, (i+1)/n]
	double c0 = grid->p[i];  // Constant coefficient in subinterval spline
	double c1 = grid->m[i];  // Coefficient of t
	double c2 = -3*c0 - 2*c1 + 3*grid->p[i+1] - grid->m[i+1]; // Coefficient of t^2
	double c3 = 2*c0 + c1 - 2*grid->p[i+1] + grid->m[i+1];    // Coefficient of t^3
	return (((c3*t+c2)*t+c1)*t+c0);
}

static inline double rbeta_4_1()
{
  int i;
  double U[4];
  for (i=0; i<4; i++)
    U[i] = rng_rand();
  for (i=0; i<3; i++) {
    if(U[i] > U[i+1]) {
      U[i+1] = U[i];
    }
  }
  return U[3]; // Largest U(0,1) of four
}

static inline double rbeta_3_2()
{
  int i, last;
  double U[4];
  for (i=0; i<4; i++)
    U[i] = rng_rand();
  for (last = 3; last >= 2; last--) {
    for (i=0; i<last; i++) {
      if(U[i] > U[i+1]) {
        double temp = U[i];
        U[i] = U[i+1];
        U[i+1] = temp;
      }
    }
  }
  return U[2]; // 2nd largest U(0,1) of four
}

static inline double rbeta_1_4() {return 1-rbeta_4_1();}
static inline double rbeta_2_3() {return 1-rbeta_3_2();}

static double spline_draw(Grid *grid)
{
	int k, n_bins = grid->n_bins, n_knots = grid->n_knots;
  double t;
	double *m = grid->m;
	double *p = grid->p;
	draw_discrete(n_knots, p, 1, &k);
	if (k==0) {
		if (rng_rand() < 3*p[0]/(6*p[0]+m[0]))
			t = rbeta_1_4();
		else
			t = rbeta_2_3();
	}
	else if (k==n_bins) {
		if (rng_rand() < 3*p[n_bins]/(6*p[n_bins]-m[n_bins]))
			t = rbeta_4_1();
		else
			t = rbeta_3_2();
	}
	else {
		t = rng_rand();
		if (rng_rand() < t*t*(3-2*t))
			t = 1-t;
		if (rng_rand() > (p[k]*(1+2*t) + m[k]*t) / (2*p[k]*(1+2*t))) {
			k = k-1;
			t = 1-t;
		}
	}
	return (k+t)/n_bins;
}

// Compute r'(x), derivative of r(x) = c exp(phi_e(x)) cosh phi_o(x)
static void compute_r_prime(int k, double *a, Grid *grid)
{
	int i;
	grid->phi_plus_prime[k] = grid->phi_minus_prime[k] = 0.0;
	for (i=1; i<=4; i++) {
		grid->phi_plus_prime[k] += a[i+1] * grid->z_plus[k + i*grid->n_knots];
		grid->phi_minus_prime[k] += a[i+1] * grid->z_minus[k + i*grid->n_knots];
	}
	double term = exp(grid->phi_plus[k]) * grid->phi_plus_prime[k];
	term -= exp(grid->phi_minus[k]) * grid->phi_minus_prime[k];
	grid->f_e_plus_prime[k] = 0.5 * term;
}

void skew_draw_eval(Skew_parameters *skew, Grid *grid)
{
  int i, k, knot, n_knots = grid->n_knots;
  double u, v, x, t;
  double phi_o; //WJM Value needs to be calculated
  double log_density;

	// Compute powers of sigma, the prior standard deviation.
	sigma[2] = skew->s2_prior;
	sigma[1] = sqrt(sigma[2]);
	sigma[3] = sigma[1] * sigma[2];
	sigma[4] = sigma[2] * sigma[2];
	sigma[5] = sigma[3] * sigma[2];

	// Compute prior-precision-normalized coefficients
	a[2] = sigma[2] * skew->h2;
	a[3] = sigma[3] * skew->h3;
	a[4] = sigma[4] * skew->h4;
	a[5] = sigma[5] * skew->h5;

	// Evalute f_even at all grid points
	grid->f_u[0] = f_u_0;
	for (knot=0; knot<=n_knots; knot++) {
		grid->phi_plus[knot] = grid->phi_minus[knot] = 0.0;
		for (i=2; i<=5; i++) {
			grid->phi_plus[knot] += a[i] * grid->z_plus[knot + i*n_knots];
			grid->phi_minus[knot] += a[i] * grid->z_minus[knot + i*n_knots];
		}
		//WJM ratio undeclared: ratio[knot] = 0.5 * (exp(grid->phi_plus[knot]) + exp(grid->phi_minus[knot]));
	}

	// Compute knot probabilities and normalization constant.
	double sum = 0.5 * f_u_0 - fp_u_0/12.0; // Contribution of first knot
	for (knot=1; knot<n_knots; knot++)
		sum += 0.0; //WJM no member named f_even in Grid grid->f_even[knot];
	sum += 0.0; // Contribution of last knot.

	// Normalize and draw or compute knot index.
	if (skew->is_draw) {
		// WJM: draw k
	}
	else {
		x = (skew->z - skew->mode)/sigma[1];
		v = 2*Phi(x)-1;
		u = F_v(v);
		k = floor(n_knots*u);
		t = n_knots*u - k;
	}

	// Evaluate derivative of f_even at this knot.
	compute_r_prime(k, a, grid);

	// Draw u, generate v then z
	if (skew->is_draw) {
		u = spline_draw(grid);
		v = inverse_F_v(u);
		x = inverse_Phi(0.5 + 0.5*v);
		if (rng_rand() * (1+exp(2*phi_o)) < 1.0)
			x = -x;
		skew->z = x*sigma[1] + skew->mode;
	}

	// Compute evaluations
	skew->log_density = log(spline_eval(u, grid));
	skew->log_density += ln_f_v(v);
	skew->log_density += log_root_2_by_pi - 0.5*x*x;
}
