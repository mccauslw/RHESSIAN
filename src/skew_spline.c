#include <stdlib.h>
#include <math.h>

#include "alias.h"
#include "skew_spline.h"
#include "skew_grid.h"
#include "RNG.h"
#include "Phi.h"
#include "spline.h"

#define max_h 5
#define max_knots 22

// Used options(digits=22), log(sqrt(2/pi)) in R
const double log_2_pi = 1.837877066409345339082;
const double p_0 = 2.0;
const double m_0K = -2.0; // Gives m_0 when m_0K/K normalization is computed.
const int n_powers = 5;
const double sign_vector[6] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
extern Grid *grids[];
extern int min_grid_points, max_grid_points;

static int K, Kp1;

// These quantities depend on the current target parameters, so cannot be
// precomputed. They have values at all grid points.
// ----------------------------------------------------------------------

// Level p of spline density at grid points and component probabilities of
// spline density
static double p[max_knots];   // Spline values (levels)
static double pi[max_knots];  // Component probabilities

// If phi(z) is a function approximating log f(y|z), then
// phi+(z) = phi(z) for z>=0, and phi-(z) = phi(-z) for z>=0.
// The next two vectors are values of these functions on a grid of (positive)
// values of z.
static double phi_plus[max_knots];  // a_2 x^2/2 + a_3 x^3/6 + a_4 x^4/24 + a_5 x^5/120
static double phi_minus[max_knots]; // a_2 x^2/2 - a_3 x^3/6 + a_4 x^4/24 - a_5 x^5/120

// Other evaluations
static double exp_phi_plus[max_knots];
static double exp_phi_minus[max_knots];
static double r[max_knots];

// Functions related to F_v(v) = 1-(1-v)^delta, for delta=0.5
static inline double inverse_F_v(double u) {return 1.0-(1.0-u)*(1.0-u);}
static inline double F_v(double v) {return 1.0 - sqrt(1.0-v);}
static inline double ln_f_v(double v) {return log(0.5) - 0.5*log(1.0-v);}

#ifndef max
static inline double max(double a, double b) {return (a>b) ? a : b;}
static inline double min(double a, double b) {return (a>b) ? b : a;}
#endif

#ifndef true
#define true (1)
#define false (0)
#endif

static double sigma[max_h+1], a[max_h+1];
static double m_0, m_K, m_Km1;

// Compute m_k. Needs to be called only for knot K-1, knots k and k+1
// satifying x_k < |x| < x_{k+1}
static inline double compute_m_k(int k, Grid *g, int is_v, int is_trunc)
{
  int i; // Index of powers of z

  // For k==0, m[k] is a known value. Quickly return.
  if (k==0) return is_v ? m_0 : 0.0;

  // Select x_plus or xu_plus grid depending on whether or not there is a
  // u = F_(v) transformation (is_v is true) or not (is_v is false).
  double *x_plus = is_v ? g->x_plus : g->xu_plus;
  double *x_minus = is_v ? g->x_minus : g->xu_minus;

  // Compute derivatives phi_plus_prime/phi_minus_prime of phi_plus/phi_minus
  double phi_plus_prime = 0.0, phi_minus_prime = 0.0;
  for (i=1; i<=4; i++) {
    phi_plus_prime += a[i+1] * x_plus[k + i*Kp1];
    phi_minus_prime -= a[i+1] * x_minus[k + i*Kp1];
  }
  double r_prime = 0.5*(exp_phi_plus[k] * phi_plus_prime
                        +exp_phi_minus[k] * phi_minus_prime);

  // Compute m_k
  double m_k;
  if (is_v) {
    m_k = r_prime * g->c[k] - p[k] * g->f_v_prime[k];
    m_k *= g->Kf_v2_inv[k];
  }
  else {
    m_k = r_prime * g->cu[k];
    m_k *= g->K_inv;
  }

  // |m_k| > 3.0 p_k is infeasible in spline draw, return feasible value if
  // is_trunc is true.
  if (is_trunc && (fabs(m_k)/p[k] > 3.0))
    m_k /= fabs(m_k)/p[k];
  return m_k;
}

// Compute odd part of phi function
static inline double phi_odd(double x, double *a, double x_bar)
{
  if (fabs(x) > x_bar) {
    return x * fabs(x) * ((a[5]/120) * x_bar * x_bar + (a[3]/6)) * x_bar;
  }
  else {
    double x2 = x*x;
    double x3 = x2 * x;
    return ((a[5]/120 * x2) + a[3]/6) * x3;
  }
}

/* If is_draw is true, draw from the skew_draw distribution and evaluate
 * the log density of the draw. If not, evaluate at the point *z
 *
 * Inputs:
 *  is_draw:       whether to draw or not
 *  n_grid_points: number of grid points in precomputed grid
 *  is_v:          whether or not to do v = F_v^{-1}(u) transformation
 *  mode:          mode of the skew_draw distribution
 *  h:             vector of derivatives (h[2] through h[5] used) of phi at mode
 *  mu:            prior mean
 *  omega:         prior precision
 *
 * Input or output according to value of id_draw:
 *  z:             pointer to value of draw/point of evaluation
 *
 * Outputs:
 *  ln_f:          pointer to log density evaluation
 */
void skew_spline_draw_eval(int is_draw, int n_grid_points, int is_v,
                           double mode, double *h, double mu, double omega,
                           double *z, double *ln_f)
{
  // Precomputed values on specified grid
  Grid *g = grids[n_grid_points];

  int i, k;           // Using i for powers, k for knots
  double x, v, u, t;  // x is z/sigma = z*sqrt(omega)

  K = g->K;
  Kp1 = K+1;
  int x_negative = 0;
  m_0 = is_v ? (m_0K * g->K_inv) : 0.0;

  // Compute powers of sigma, the prior standard deviation
  sigma[2] = 1.0/omega;
  sigma[1] = sqrt(sigma[2]);
  sigma[3] = sigma[1] * sigma[2];
  sigma[4] = sigma[2] * sigma[2];
  sigma[5] = sigma[3] * sigma[2];

  // Compute prior-precision-normalized coefficients a[2] through a[5]
  // Note that these give the same value:
  //   a[2] x^2/2 + a[3] x^3/6 + a[4] x^4/24 + a[5] x^5/120
  //   h[2] z^2/2 + h[3] z^3/6 + h[4] z^4/24 + h[5] z^5/120
  for (i=2; i<=5; i++)
    a[i] = sigma[i] * h[i];
  double x_bar = exp(log(12.0/fabs(a[3]))/3.0);

  // Evalute level p[k] = f_even at all grid points k = 0, 1, ..., K
  double *x_plus = is_v ? g->x_plus : g->xu_plus;
  for (k=1; k<=K; k++) {
    double plus_term = a[2] * x_plus[k + 2*Kp1] + a[4] * x_plus[k + 4*Kp1];
    double alt_term = a[3] * x_plus[k + 3*Kp1] + a[5] * x_plus[k + 5*Kp1];
    phi_plus[k] = plus_term + alt_term;
    phi_minus[k] = plus_term - alt_term;
    exp_phi_plus[k] = exp(phi_plus[k]);
    exp_phi_minus[k] = exp(phi_minus[k]);
    p[k] = r[k] = 0.5 * (exp_phi_plus[k] + exp_phi_minus[k]);
    if (is_v) p[k] *= g->f_v_inv[k];
  }

  // Computation of p[K] and m[K] are based on true computed values
  double c_p = 0.1, c_m = 0.1;
  double p_tau = p[K], m_tau = compute_m_k(K, g, is_v, false);
  m_Km1 = compute_m_k(K-1, g, is_v, false);
  if (is_v) {
    p[0] = p_0;
    p[K] = 0.0;
    m_K = min(c_m*m_tau, -6.0*p[K-1] - 4.0*m_tau - m_Km1);
  }
  else {
    p[0] = 1.0;
    p[K] = -4.0*p_tau + 5.0*p[K-1] + 2*m_tau + m_Km1;
    m_K = min(0.0, 24*(p[K-1]-p_tau) + 8.0*m_tau + 5.0*m_Km1);
    if (p[K] < 0.0) {
      double p_min = c_p*p_tau;
      m_K += 6.0*(p_min - p[K]);
      p[K] = p_min;
    }
    m_K = min(0.0, m_K);
  }

  // From here on we need to truncate to avoid negative density values.
  if (m_Km1 < -3.0*p[K-1])
    m_Km1 = -3.0*p[K-1];

  // Compute knot probabilities and normalization constant.
  double p_Delta = p_tau - 0.5*(p[K-1] + p[K]) - 0.125*(m_Km1 - m_K);
  double m_Delta = m_tau - 1.5*(p[K] - p[K-1]) + 0.25*(m_Km1 + m_K);
  if (p_Delta < 0.0) {
    p_Delta = 0.0;
    m_Delta = 0.0;
  }
  else if (fabs(m_Delta) > 3.0*p_Delta) {
    m_Delta *= 3.0*p_Delta/fabs(m_Delta);
  }

  double pi_total = (pi[0] = 0.5*p[0] + m_0/12); // Contribution of first knot
  for (k=1; k<K; k++)
    pi_total += (pi[k] = p[k]);
  pi_total += (pi[K] = 0.5*p[K] - m_K/12);       // Contribution of last knot
  pi_total += (pi[K+1] = 0.5*p_Delta);           // Contribution of tau knot

  static int done = 0;
  if (!done) {
    for (k=0; k<=K+1; k++)
      printf("k=%d, pi[k]=%lf, p[k]=%lf\n", k, pi[k], p[k]);
    printf("pi_total=%lf\n", pi_total);

    printf("p[K-1] = %le, m[K-1] = %le\np_tau = %le, m_tau = %le\n"
             "p[K] = %le, m[K] = %le\np_Delta = %le, m_Delta = %le\n",
           p[K-1], m_Km1, p_tau, m_tau, p[K], m_K, p_Delta, m_Delta);
    printf("a2 = %lf, a3 = %lf, a4 = %lf, a5 = %lf\n", a[2], a[3], a[4], a[5]);
    printf("x_bar = %lf\n", x_bar);
    done = 1;
  }

  // Repeatable part of draw, in case a loop is desired.
  double phi_o, exp_2_phi_o; // Odd part of phi function and exp of twice odd part

  // Normalize and draw or compute knot index.
  if (is_draw)
    draw_discrete(K+2, pi, 1, &k);
  else {
    x = (*z - mode)/sigma[1];
    phi_o = phi_odd(x, a, x_bar);
    exp_2_phi_o = exp(2.0 * phi_o);
    x_negative = (x < 0);
    x = fabs(x);
    v = (x==0) ? 0 : 2*Phi(x)-1;
    u = is_v ? F_v(v) : v;
    k = floor(K*u);
    t = K*u - k;
  }

  // Evaluate derivative of f_even at this knot.
  double m_k, m_kp1;
  if (k<K-1)
    m_k = compute_m_k(k, g, is_v, true);
  else
    m_k = (k==K-1) ? m_Km1 : m_K;

  // Draw u, generate v then z
  if (is_draw) {

    // Draw u
    if (k==0) {
      t = left_t_draw(p[0], m_k);
      m_kp1 = compute_m_k(1, g, is_v, true);
    }
    else if (k==K) {
      k = k-1;
      m_k = m_Km1;
      m_kp1 = m_K;
      t = right_t_draw(p[K], m_K);
    }
    else if (k==K+1) {
      k = k-2;
      m_k = m_Km1;
      m_kp1 = m_K;
      t = (inner_t_draw(p_Delta, m_Delta/2) + 1) / 2;
    }
    else {
      t = inner_t_draw(p[k], m_k);
      if (t<0) {
        t = t+1; k = k-1;
        m_kp1 = m_k;
        m_k = compute_m_k(k, g, is_v, true);
      }
      else {
        m_kp1 = (k==K-1) ? m_K : compute_m_k(k+1, g, is_v, true);
      }
    }
    u = (k+t)/K;
    v = is_v ? inverse_F_v(u) : u;
    x = inverse_Phi(0.5 + 0.5*v);
    phi_o = phi_odd(x, a, x_bar);
    if (rng_rand() * cosh(phi_o) * 2 < exp(-phi_o)) {
      x = -x;
      phi_o = -phi_o;
    }
    *z = x*sigma[1] + mode;
  }
  else {
    m_kp1 = (k==K-1) ? m_K : ((k==K-2) ? m_Km1 : compute_m_k(k+1, g, is_v, true));
  }

  // Compute exact value of log approximate density
  // ----------------------------------------------

  // Compute f_u factor
  double f_u, c0, c1, c2, c3;
  c0 = p[k];          // Constant coefficient in subinterval spline
  c1 = m_k;           // Coefficient of t
  c2 = -3*c0 - 2*c1 + 3*p[k+1] - m_kp1; // Coefficient of t^2
  c3 = 2*c0 + c1 - 2*p[k+1] + m_kp1;    // Coefficient of t^3
  f_u = (((c3*t+c2)*t+c1)*t+c0);
  if (k==(K-1)) {
    double t_tilde = 2*t-1; // In [-1, 1], corresponding to [u[K-1], 1]
    double tt_abs = fabs(t_tilde);
    double tt_sign = (t_tilde > 0) ? 1.0 : -1.0;
    f_u += ((2*tt_abs - 3)*tt_abs*tt_abs + 1) * p_Delta;
    f_u += tt_sign * (((tt_abs - 2.0)*tt_abs + 1.0)*tt_abs) * 0.5 * m_Delta;
  }

  // Compute terms of log approximate density
  *ln_f = log(f_u) - log(pi_total) + g->log_K;
  if (is_v)
    *ln_f += ln_f_v(v);
  *ln_f += 0.5 * (log(omega) - log_2_pi - x*x);
  *ln_f += phi_o - log(cosh(phi_o));
}
