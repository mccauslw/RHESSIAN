#include <R.h>
#include <Rinternals.h>

#include <stdlib.h>
#include <math.h>
#include "new_skew.h"
#include "alias.h"
#include "spline.h"
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
		spline_draw(grid->n_knots, grid->p, grid->m, 1, &u);
		v = inverse_F_v(u);
		x = inverse_Phi(0.5 + 0.5*v);
		if (rng_rand() * (1+exp(2*phi_o)) < 1.0)
			x = -x;
		skew->z = x*sigma[1] + skew->mode;
	}

	// Compute evaluations
	double f_u;
	spline_eval(grid->n_knots, grid->p, grid->m, 1, &u, &f_u);
	skew->log_density = log(f_u);
	skew->log_density += ln_f_v(v);
	skew->log_density += log_root_2_by_pi - 0.5*x*x;
}
