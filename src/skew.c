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
const double p_0 = 2.0;
const double m_0K = -2.0; // Gives m_0 when m_0K/K normalization is computed.
const int n_powers = 5;
extern Grid g6;

/* These depend on the current target parameters, and have values at all grid
 * points.
 */
// p and m and parameters of the spline density
static double p[max_knots]; // Spline values (levels)
static double m[max_knots]; // Spline derivatives, normalized for t in [0,1]
static double pi[max_knots]; // Component probabilities

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

// Functions related to F_v(v) = 1-(1-v)^delta, for delta=0.5
static inline double inverse_F_v(double u) {return 1.0-(1.0-u)*(1.0-u);}
static inline double F_v(double v) {return 1.0 - sqrt((1.0-v));}
static inline double f_v(double v) {return 0.5 / sqrt(1.0-v);}
static inline double ln_f_v(double v) {return log(0.5) - 0.5*log(1.0-v);}

#ifndef max
static inline double max(double a, double b) {return (a>b) ? a : b;}
#endif

static double sigma[max_h+1], a[max_h+1];

// Compute m_k. Needs to be called only for second last knot and value that is
// actually drawn.
static double compute_m_k(int k, Grid *g)
{
  int K = g->K;
  int i; // Index of powers of z
  // phi'(z) = phi'+(z) when z is postive, phi'(z) = phi'-(-z) when z is negative
  // phi_plus_prime and phi_minus_prime are the derivatives of phi_plus, phi_minus
  double phi_plus_prime = 0.0, phi_minus_prime = 0.0;
  if (k==0) return m_0K / g->K;
  for (i=1; i<=4; i++) {
    phi_plus_prime += a[i+1] * g->x_plus[k + i*(K+1)];
    phi_minus_prime -= a[i+1] * g->x_minus[k + i*(K+1)];
  }
  double r_prime = 0.5*(exp_phi_plus[k] * phi_plus_prime
                        +exp_phi_minus[k] * phi_minus_prime);
  double m_k = r_prime * g->c[k] - p[k] * g->f_v_prime[k];
  m_k /= g->K * g->f_v[k] * g->f_v[k];
  return m_k;
}

// Compute odd part of phi function
static inline double phi_odd(double z, double *a)
{
  double z2 = z*z;
  double z3 = z2 * z;
  return ((a[5]/120 * z2) + a[3]/6) * z3;
}

/* If is_draw is true, draw from the skew_draw distribution and evaluate
 * the log density of the draw. If not, evaluate at the point *z
 *
 * Inputs:
 *  is_draw:  whether to draw or not
 *  n_draws:  number of draws/points of evaluation
 *  g:        precomputed grid to use
 *  mode:     mode of the skew_draw distribution
 *  h:        vector of derivatives (h[2] through h[5] used) of phi at mode
 *  mu:       prior mean
 *  omega:    prior precision
 *
 * Input or output according to value of id_draw:
 *  z:        vector of values
 *
 * Outputs:
 *  ln_f:     vector of log density evaluations
 */
void skew_draw_eval(int is_draw, int n_draws, Grid *g,
                    double mode, double *h, double mu, double omega,
                    double *z, double *ln_f)
{
  int i, k, draw, K = g->K; // Using i for powers, k for knots
  double x, v, u, t;
  int x_negative = 0;

  // Compute powers of sigma, the prior standard deviation
  sigma[2] = 1.0/omega;
  sigma[1] = sqrt(sigma[2]);
  sigma[3] = sigma[1] * sigma[2];
  sigma[4] = sigma[2] * sigma[2];
  sigma[5] = sigma[3] * sigma[2];

  // Compute prior-precision-normalized coefficients
  for (i=2; i<=5; i++)
    a[i] = sigma[i] * h[i];

  // Evalute level p[k] = f_even at all grid points k=0,1,...,K
  p[0] = p_0;
  for (k=1; k<=K; k++) {
    phi_plus[k] = phi_minus[k] = 0.0;
    for (i=2; i<=5; i++) {
      phi_plus[k] += a[i] * g->x_plus[k + i*(K+1)];
      phi_minus[k] += a[i] * g->x_minus[k + i*(K+1)];
      exp_phi_plus[k] = exp(phi_plus[k]);
      exp_phi_minus[k] = exp(phi_minus[k]);
    }
    r[k] = 0.5 * (exp_phi_plus[k] + exp_phi_minus[k]);
    p[k] = r[k] / g->f_v[k];
  }

  // Compute knot probabilities and normalization constant.
  double m_Km1 = compute_m_k(K-1, g);
  double p_Delta = max(0, p[K] - 0.5*p[K-1] - 0.125*m_Km1);
  double pi_total = (pi[0] = p[0]/2 - m[0]/12); // Contribution of first knot
  for (k=1; k<K; k++)
    pi_total += (pi[k] = p[k]);
  pi_total += (pi[K] = p_Delta/2); // Contribution of last knot.

  for (draw=0; draw<n_draws; draw++) {
    double phi_o; // Odd part of phi function

    // Normalize and draw or compute knot index.
    if (is_draw)
      draw_discrete(K+1, pi, 1, &k);
    else {
      x = (z[draw] - mode)/sigma[1];
      phi_o = phi_odd(x, a);
      x_negative = (x < 0);
      x = fabs(x);
      v = (x==0) ? 0 : 2*Phi(x)-1;
      u = (x==0) ? 0 : F_v(v);
      k = floor(K*u);
      t = K*u - k;
    }

    // Evaluate derivative of f_even at this knot.
    double m_k = compute_m_k(k, g);
    double m_kp1 = compute_m_k(k+1, g);

    // Draw u, generate v then z
    if (is_draw) {

      // Draw u
      if (k==0)
        t = left_t_draw(p[0], m[0]);
      else if (k==K) {
        double m_Delta = (m_k + 1.5 * p[K-1] + 0.25 * m_Km1) / 2;
        k = k-1;
        t = (inner_t_draw(p_Delta, m_Delta) + 1) / 2;
      }
      else {
        t = inner_t_draw(p[k], m[k]);
        if (t<0) {
          t = t+1; k = k-1;
          m_kp1 = m_k;
          m_k = compute_m_k(k, g);
        }
      }
      u = (k+t)/K;
      v = inverse_F_v(u);
      x = inverse_Phi(0.5 + 0.5*v);
      phi_o = phi_odd(x, a);
      if (rng_rand() * (1+exp(2*phi_o)) < 1.0)
        x = -x;
      z[draw] = x*sigma[1] + mode;
    }

    // Compute evaluations
    double c0, c1, c2, c3, f_u;
    c0 = p[k];          // Constant coefficient in subinterval spline
    c1 = m_k;           // Coefficient of t
    c2 = -3*c0 - 2*c1 + 3*p[k+1] - m_kp1; // Coefficient of t^2
    c3 = 2*c0 + c1 - 2*p[k+1] + m_kp1;    // Coefficient of t^3
    f_u = (((c3*t+c2)*t+c1)*t+c0);

    //spline_eval(K, p, m, 1, &u, &f_u);
    ln_f[draw] = log(f_u) - log(pi_total);
    ln_f[draw] += ln_f_v(v);
    ln_f[draw] += log(1.0 + (exp(2*phi_o)-1) / (exp(2*phi_o)+1));
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
  skew_draw_eval(FALSE, length(z), &g6, mode, h_ptr, mu, omega, z_ptr, log_f_ptr);
  UNPROTECT(3);
  return log_f;
}

SEXP skew_draw_c(double mode, SEXP h, double mu, double omega, int n_draws)
{
  double *h_ptr, *draws_ptr, *log_f_ptr;
  SEXP draws, log_f;
  h = PROTECT(coerceVector(h, REALSXP));
  draws = PROTECT(allocVector(REALSXP, n_draws));
  log_f = PROTECT(allocVector(REALSXP, n_draws));
  h_ptr = REAL(h);
  draws_ptr = REAL(draws);
  log_f_ptr = REAL(log_f);
  skew_draw_eval(TRUE, n_draws, &g6, mode, h_ptr, mu, omega, draws_ptr, log_f_ptr);
  UNPROTECT(3);
  return draws;
}
