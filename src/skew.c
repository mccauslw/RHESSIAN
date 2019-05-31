#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>

#include "alias.h"
#include "skew.h"
#include "skew_grid.h"
#include "RNG.h"
#include "Phi.h"
#include "spline.h"

#define max_h 5
#define max_knots 20

// Used options(digits=22), log(sqrt(2/pi)) in R
const double log_root_2_by_pi = -0.2257913526447273833142;
const double f_u_0 = 2.0;
const double fp_u_0 = -2.0;
const int n_powers = 5;
extern Grid g6;

/* These depend on the current target parameters, and have values at all grid
 * points.
 */
// p and m and parameters of the spline density
static double p[max_knots]; // Spline values (levels)
static double m[max_knots]; // Spline derivatives, normalized for t in [0,1]

// phi(z) is a function approximating the log of the target density.
// phi(z) = phi+(z) when z is positive, phi-(-z) when z is negative
// The next two vectors are values of these functions on a grid of (positive)
// values of z.
static double phi_plus[max_knots];  // a_2 z^2/2 + a_3 z^3/6 + a_4 z^4/24 + a_5 z^5/120 on grid
static double phi_minus[max_knots]; // a_2 z^2/2 - a_3 z^3/6 + a_4 z^4/24 - a_5 z^5/120 on grid

// Other evaluations
static double exp_phi_plus[max_knots];
static double exp_phi_minus[max_knots];
static double r[max_knots];
static double r_prime[max_knots];

//static double f_e_plus_prime[max_knots];
//static double f_u[max_knots];
//static double ratio[max_knots];
//static double f_even[max_knots];

// Functions related to F_v(v) = 1-(1-v)^delta, for delta=0.5
static inline double inverse_F_v(double u) {return 1.0-(1.0-u)*(1.0-u);}
static inline double F_v(double v) {return 1.0 - sqrt((1.0-v));}
static inline double f_v(double v) {return 0.5 / sqrt(1.0-v);}
static inline double ln_f_v(double v) {return log(0.5) - 0.5*log(1.0-v);}

static void grid_initialize(int use_delta, Grid *g)
{
  int K = g->K;
  int k, i;
  for (k=0; k<K; k++) {
    // Compute u and various transformation of it at the grid points
    g->u[k] = (k==K-1) ? ((double) k)/K : 1 - 0.5/K;
    g->v[k] = inverse_F_v(g->u[k]);
    g->z[k] = inverse_Phi(0.5 + 0.5*g->v[k]);
    g->f_v[k] = f_v(g->v[k]);
    g->f_v_prime[k] = 0.5 * g->f_v[k] / (1-g->v[k]);

    // Precompute powers of z and -z at current knot.
    g->z_plus[k] = g->z_minus[k] = 1.0;
    for (i=1; i<=n_powers; i++) {
      g->z_plus[k + i*K] = g->z_plus[k + (i-1)*K] * g->z[k] / i;
      g->z_minus[k + i*K] = g->z_minus[k + (i-1)*K] * -g->z[k] / i;
    }
  }
}

static double sigma[max_h+1], a[max_h+1];

// Compute r'(x), derivative of r(x) = c exp(phi_e(x)) cosh phi_o(x)
// This function needs to be called only for the knot value that is
// actually drawn.
static double compute_m_k(int k, Grid *g)
{
  int K = g->K;
  int i; // Index of powers of z
  // phi'(z) = phi'+(z) when z is postive, phi'(z) = phi'-(-z) when z is negative
  // phi_plus_prime and phi_minus_prime are the derivatives of phi_plus, phi_minus
  double phi_plus_prime = 0.0, phi_minus_prime = 0.0;
  for (i=1; i<=4; i++) {
    phi_plus_prime += a[i+1] * g->z_plus[k + i*K];
    phi_minus_prime += a[i+1] * g->z_minus[k + i*K];
  }
  double r_prime = 0.5*(exp_phi_plus[k] * phi_plus_prime
                        -exp_phi_minus[k] * phi_minus_prime);
  double m_k = r_prime - p[k] * g->f_v_prime[k];
  m_k /= g->K * g->f_v[k] * g->f_v[k];
  return m_k;
}

void skew_draw_eval(double mode, double *h, double mu, double omega,
                    int n, double *z, double *ln_f,
                    int is_draw, Grid *g)
{
  int i, k, draw, K = g->K; // Using i for powers, k for knots
  double x, v, u, t;
  grid_initialize(1, g);

  // Compute powers of sigma, the prior standard deviation
  sigma[2] = 1.0/omega;
  sigma[1] = sqrt(sigma[2]);
  sigma[3] = sigma[1] * sigma[2];
  sigma[4] = sigma[2] * sigma[2];
  sigma[5] = sigma[3] * sigma[2];

  // Compute prior-precision-normalized coefficients
  for (i=2; i<=5; i++)
    a[i] = sigma[i] * h[i];

  // Evalute level p[k] f_even at all grid points k=0,1,...,K
  p[0] = f_u_0;
  for (k=0; k<=K; k++) {
    phi_plus[k] = phi_minus[k] = 0.0;
    for (i=2; i<=5; i++) {
      phi_plus[k] += a[i] * g->z_plus[k + i*K];
      phi_minus[k] += a[i] * g->z_minus[k + i*K];
      exp_phi_plus[k] = exp(phi_plus[k]);
      exp_phi_minus[k] = exp(phi_minus[k]);
    }
    r[k] = 0.5 * (exp_phi_plus[k] + exp_phi_minus[k]);
    p[k] = r[k] / g->f_v[k];
  }

  // Compute knot probabilities and normalization constant.
  double sum = 0.5 * f_u_0 - fp_u_0/12.0; // Contribution of first knot
  for (k=1; k<K; k++)
    sum += p[k];
  sum += 0.0; // Contribution of last knot.

  for (draw=0; draw<n; draw++) {
    double phi_o;

    // Normalize and draw or compute knot index.
    if (is_draw) {
      // WJM: draw k
    }
    else {
      x = (z[draw] - mode)/sigma[1];
      v = 2*Phi(x)-1;
      u = F_v(v);
      k = floor(K*u);
      t = K*u - k;
    }

    // Evaluate derivative of f_even at this knot.
    double m_k = compute_m_k(k, g);

    // Draw u, generate v then z
    if (is_draw) {
      spline_draw(K, p, m, 1, &u);
      v = inverse_F_v(u);
      x = inverse_Phi(0.5 + 0.5*v);
      if (rng_rand() * (1+exp(2*phi_o)) < 1.0)
        x = -x;
      z[draw] = x*sigma[1] + mode;
    }

    // Compute evaluations
    double f_u;
    spline_eval(K, p, m, 1, &u, &f_u);
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
  skew_draw_eval(mode, h_ptr, mu, omega, length(z), z_ptr, log_f_ptr, FALSE, &g6);
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
  skew_draw_eval(mode, h_ptr, mu, omega, n_draws, draws_ptr, NULL, TRUE, &g6);
  UNPROTECT(2);
  return draws;
}
