#ifndef SKEW
#define SKEW

/*
 * The Skew_parameters structure ...
 */
typedef struct {
  double mode;
  double h2;
  double h3;
  double h4;
  double h5;
  double s2_prior;
  double u_sign;
  int is_draw;
  double z;
  int n_reject;
  double log_density;
} Skew_parameters;

/*
  The Grid structure gathers together the evaluation of several functions on a grid
  of points (u_0,...,u_n) = (0, 1/n, 2/n, ..., (n-1)/n, 1-1/(2n)). Note that except
  for the last one, the grid points are uniform.
  The functions are related to a cubic spline approximation of a density on [0, 1].

  See documentation in "Flexible even densities", file name flexeven.Rmd for details.
*/
typedef struct {
	int n_knots;       // Number of spline knots
	int n_bins;        // Number of spline bins
	int use_delta;     // Whether (1) or not (0) to use delta=0.5 transformation

	double *u;         // Values of u on grid
	double *v;         // Values of v on grid, where v = F_v^{-1}(u)
	double *z;         // Values of z on grid, where z = Phi^{-1}((1+v/2))

	double *f_v;       // Density f_v of v on grid, where F_v is cdf of v
	double *f_v_prime; // Derivative of density f_v of v on grid

	double *p;         // Spline values (levels)
	double *m;         // Spline derivatives, normalized for t in [0,1]

	// phi_plus and phi_minus are even and odd, respectively, parts of the function
	// phi approximating the log of the target density f() we want to approximate.
	double *phi_plus;  // a_2 z^2/2 + a_3 z^3/6 + a_4 z^4/24 + a_5 z^5/120 on grid
	double *phi_minus; // a_2 z^2/2 - a_3 z^3/6 + a_4 z^4/24 - a_5 z^5/120 on grid

	// phi_plus_prime and phi_minus_prime are the derivatives of phi_plus, phi_minus
	double *phi_plus_prime;  //  a_2 z + a_3 z^2/2 + a_4 z^3/6 + a_5 z^4/24 on grid
	double *phi_minus_prime; // -a_2 z + a_3 z^2/2 - a_4 z^3/6 + a_5 z^4/24 on grid

	double *f_e_plus;        // Positive density of even part of target
	double *f_e_plus_prime;  // Derivative of f_e_plus

	double *z_plus;          // Powers i=0,1,2,3,4,5 of z_k, k=0,1,...,n
	double *z_minus;         // Powers i=0,1,2,3,4,5 of -z_k, k=0,1,...,n

	// Might be redundant
	double *f_u;
	double *fp_u;
} Grid;

#endif
