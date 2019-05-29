#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>

#include "alias.h"
#include "skew.h"
#include "RNG.h"
#include "Phi.h"
#include "spline.h"

#define max_h 5
#define max_knots 20

static int global_n_knots = 0;
static int use_delta = 1; // True

// Used options(digits=22), log(sqrt(2/pi)) in R
const double log_root_2_by_pi = -0.2257913526447273833142;
const double delta = 0.5;
const double delta_dual = 2.0;
const double f_u_0 = delta_dual;
const double fp_u_0 = -delta * delta_dual * delta_dual;

// These depend only on n_knots and use_delta
static double u[max_knots];
static double v[max_knots];
static double z[max_knots];
static double f_v[max_knots];
static double f_v_prime[max_knots];
static double z_plus[max_h+1];
static double z_minus[max_h+1];

/* These depend on the current target parameters, and have values at all grid points
 */
// p and m and parameters of the spline density
static double p[max_knots]; // Spline values (levels)
static double m[max_knots]; // Spline derivatives, normalized for t in [0,1]

// phi_plus and phi_minus are even and odd, respectively, parts of the function
// phi approximating the log of the target density f() we want to approximate.
static double phi_plus[max_knots];  // a_2 z^2/2 + a_3 z^3/6 + a_4 z^4/24 + a_5 z^5/120 on grid
static double phi_minus[max_knots]; // a_2 z^2/2 - a_3 z^3/6 + a_4 z^4/24 - a_5 z^5/120 on grid

// phi_plus_prime and phi_minus_prime are the derivatives of phi_plus, phi_minus
static double phi_plus_prime[max_knots];  //  a_2 z + a_3 z^2/2 + a_4 z^3/6 + a_5 z^4/24 on grid
static double phi_minus_prime[max_knots]; // -a_2 z + a_3 z^2/2 - a_4 z^3/6 + a_5 z^4/24 on grid

// Other evaluations
static double f_e_plus_prime[max_knots]; // Even part of target density
static double f_u[max_knots]; //WJM: comment these
static double ratio[max_knots];
static double f_even[max_knots];

// Functions related to F_v(v) = 1-(1-v)^delta, for delta=0.5
static inline double inverse_F_v(double u) {return 1.0-(1.0-u)*(1.0-u);}
static inline double F_v_fcn(double v) {return 1.0-sqrt(1.0-v);}
static inline double f_v_fcn(double v) {return 0.5/sqrt(1.0-v);}
static inline double ln_f_v(double v) {return log(0.5) - 0.5*log(1.0-v);}

static void grid_initialize(int n_knots, int n_powers, int use_delta)
{
  global_n_knots = n_knots;
  int k, i;
  for (k=0; k<n_knots; k++) {
    // Compute u and various transformation of it at the grid points
    u[k] = (k==n_knots-1) ? ((double) k)/n_knots : 1 - 0.5/n_knots;
    v[k] = inverse_F_v(u[k]);
    z[k] = inverse_Phi(0.5 + 0.5*v[k]);
    f_v[k] = f_v_fcn(v[k]);
    f_v_prime[k] = 0.5*f_v[k]/(1-v[k]);

    // Precompute powers of z and -z at current knot.
    z_plus[k] = z_minus[k] = 1.0;
    for (i=1; i<=n_powers; i++) {
      z_plus[k + i*n_knots] = z_plus[k + (i-1)*n_knots] * z[k] / i;
      z_minus[k + i*n_knots] = z_minus[k + (i-1)*n_knots] * -z[k] / i;
    }
  }
}

static double sigma[max_h+1], a[max_h+1];

// Compute r'(x), derivative of r(x) = c exp(phi_e(x)) cosh phi_o(x)
static void compute_r_prime(int k, int n_knots)
{
  int i;
  phi_plus_prime[k] = phi_minus_prime[k] = 0.0;
  for (i=1; i<=4; i++) {
    phi_plus_prime[k] += a[i+1] * z_plus[k + i*n_knots];
    phi_minus_prime[k] += a[i+1] * z_minus[k + i*n_knots];
  }
  double term = exp(phi_plus[k]) * phi_plus_prime[k];
  term -= exp(phi_minus[k]) * phi_minus_prime[k];
  f_e_plus_prime[k] = 0.5 * term;
}

void skew_draw_eval(double mode, double *h, double mu, double omega,
                    int n, double *z, double *ln_f, int is_draw, int n_knots)
{
  int i, k, draw; // Using i for powers, k for knots
  if (global_n_knots == 0)
    grid_initialize(8, 5, 1);

  // Compute powers of sigma, the prior standard deviation
  sigma[2] = 1.0/omega;
  sigma[1] = sqrt(sigma[2]);
  sigma[3] = sigma[1] * sigma[2];
  sigma[4] = sigma[2] * sigma[2];
  sigma[5] = sigma[3] * sigma[2];

  // Compute prior-precision-normalized coefficients
  for (i=2; i<=5; i++)
    a[i] = sigma[i] * h[i];

  // Evalute f_even at all grid points
  f_u[0] = f_u_0;
  for (k=0; k<=n_knots; k++) {
    phi_plus[k] = phi_minus[k] = 0.0;
    for (i=2; i<=5; i++) {
      phi_plus[k] += a[i] * z_plus[k + i*n_knots];
      phi_minus[k] += a[i] * z_minus[k + i*n_knots];
    }
    ratio[k] = 0.5 * (exp(phi_plus[k]) + exp(phi_minus[k]));
  }

  // Compute knot probabilities and normalization constant.
  double sum = 0.5 * f_u_0 - fp_u_0/12.0; // Contribution of first knot
  for (k=1; k<n_knots; k++)
    sum += f_even[k];
  sum += 0.0; // Contribution of last knot.

  for (draw=0; draw<n; draw++) {
    double phi_o;

    // Normalize and draw or compute knot index.
    if (is_draw) {
      // WJM: draw k
    }
    else {
      x = (z - mode)/sigma[1];
      v = 2*Phi(x)-1;
      u = F_v(v);
      k = floor(n_knots*u);
      t = n_knots*u - k;
    }

    // Evaluate derivative of f_even at this knot.
    compute_r_prime(k, n_knots);

    // Draw u, generate v then z
    if (is_draw) {
      spline_draw(n_knots, p, m, 1, &u);
      v = inverse_F_v(u);
      x = inverse_Phi(0.5 + 0.5*v);
      if (rng_rand() * (1+exp(2*phi_o)) < 1.0)
        x = -x;
      z = x*sigma[1] + mode;
    }

    // Compute evaluations
    double f_u;
    spline_eval(n_knots, p, m, 1, &u, &f_u);
    ln_f[draw] = log(f_u);
    ln_f[draw] += ln_f_v(v);
    ln_f[draw] += log_root_2_by_pi - 0.5*x*x;
  }
}

SEXP skew_eval_c(double mode, SEXP h, double mu, double omega, SEXP z)
{
  double *h_ptr, *z_ptr, *log_f_ptr;
  SEXP log_f;
  h = PROTECT(coerceVector(h, REALSXP));
  z = PROTECT(coerceVector(z, REALSXP));
  log_f = PROTECT(allocVector(REALSXP, length(z)));
  h_ptr = REAL(h);
  z_ptr = REAL(z);
  log_f_ptr = REAL(log_f);
  skew_draw_eval(mode, h_ptr, mu, omega, length(z), z_ptr, log_f_ptr);
  UNPROTECT(3);
  return log_f;
}

SEXP skew_draw_c(double mode, SEXP h, double mu, double omega, int n_draws)
{
  double *h_ptr, *draws_ptr;
  SEXP draws;
  h = PROTECT(coerceVector(h, REALSXP));
  draws = PROTECT(allocVector(REALSXP, n_draws));
  h_ptr = REAL(h);
  draws_ptr = REAL(draws);
  skew_draw_eval(mode, h_ptr, mu, omega, n_draws, draws_ptr, NULL);
  UNPROTECT(2);
  return draws;
}
