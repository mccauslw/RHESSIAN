#include <R.h>
#include <Rinternals.h>

#include "spline.h"
#include "alias.h"
#include "RNG.h"

#define max_n 100

// Specialized beta draws for select small integral parameters

// Draw a beta(4, 1) variate
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

// Draw a beta(3, 2) variate
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

static inline double rbeta_1_4() {return 1-rbeta_4_1();} // Draw a beta(1, 4) variate
static inline double rbeta_2_3() {return 1-rbeta_3_2();} // Draw a beta(2, 3) variate

void spline_eval(int n_knots, double *p, double *m, int n_evals, double *u, double *f_u)
{
  int i, b, K = n_knots-1;
  double c0, c1, c2, c3, t;
  for (i=0; i<n_evals; i++) {
    b = floor(u[i]*K);  // Bin (subinterval) index in {0,1,...,n-1}
    t = u[i]*K - b;     // index of u in subinterval [i/n, (i+1)/n]
    c0 = p[b];          // Constant coefficient in subinterval spline
    c1 = m[b];          // Coefficient of t
    c2 = -3*c0 - 2*c1 + 3*p[b+1] - m[b+1]; // Coefficient of t^2
    c3 = 2*c0 + c1 - 2*p[b+1] + m[b+1];    // Coefficient of t^3
    f_u[i] = (((c3*t+c2)*t+c1)*t+c0);
  }
}

double inner_t_draw(double p_k, double m_k)
{
  double t = rng_rand();
  double U = rng_rand();
  double Pr_comp = t*t*(3-2*t);
  double positive_threshold;
  // This code avoids an extra U(0, 1) draw.
  // Equivalent procedure that is easier to understand is
  // (1) Replace t = 1-t with probability t*t*(3-2*t). (otherwise leave it alone).
  // (2) Keep positive with probability (p_k*(3-2*t) + m_k*(1-t)) / (2*p_k*(3-2*t))
  //     computed *after* possible replacement of t by 1-t. (otherwise negate t)
  if (U < Pr_comp) {
    t = 1-t;
    positive_threshold = Pr_comp * (p_k*(1+2*t) + m_k*t) / (2*p_k*(1+2*t));
  }
  else
    positive_threshold = Pr_comp + (1-Pr_comp) * (p_k*(1+2*t) + m_k*t) / (2*p_k*(1+2*t));
  return (U < positive_threshold) ? t : -t;
}

double left_t_draw(double p_0, double m_0)
{
  return (rng_rand() < 3*p_0/(6*p_0+m_0)) ? rbeta_1_4() : rbeta_2_3();
}

double right_t_draw(double p_K, double m_K)
{
  return (rng_rand() < 3*p_K/(6*p_K-m_K)) ? rbeta_4_1() : rbeta_3_2();
}

void spline_draw(int n_knots, double *p, double *m, int n_draws, double *u)
{
  double pmf[max_n];
  int Alias[max_n];
  double Prob[max_n];
  int i, k, K = n_knots-1;
  memcpy(pmf, p, n_knots*sizeof(double));
  pmf[0] = p[0]/2 + m[0]/12;
  pmf[K] = p[K]/2 - m[K]/12;
  alias_tables(n_knots, pmf, Alias, Prob);
  for (i=0; i<n_draws; i++) {
    double t;
    draw_discrete_from_alias_tables(n_knots, Alias, Prob, 1, &k); // Draw random knot
    if (k==0)
      t = left_t_draw(p[0], m[0]);
    else if (k==K) {
      k = k-1;
      t = right_t_draw(p[K], m[K]);
    }
    else
      t = inner_t_draw(p[k], m[k]); // t may be negative (projecting into previous bin)
    u[i] = (k+t)/K;
  }
}

SEXP spline_eval_c(SEXP p, SEXP m, SEXP u)
{
  double *p_ptr, *m_ptr, *u_ptr, *f_u_ptr;
  SEXP f_u;
  p = PROTECT(coerceVector(p, REALSXP));
  m = PROTECT(coerceVector(m, REALSXP));
  u = PROTECT(coerceVector(u, REALSXP));
  f_u = PROTECT(allocVector(REALSXP, length(u)));
  p_ptr = REAL(p);
  m_ptr = REAL(m);
  u_ptr = REAL(u);
  f_u_ptr = REAL(f_u);
  spline_eval(length(p), p_ptr, m_ptr, length(u), u_ptr, f_u_ptr);
  UNPROTECT(4);
  return f_u;
}

SEXP spline_draw_c(SEXP p, SEXP m, int n_draws)
{
  double *p_ptr, *m_ptr, *draws_ptr;
  SEXP draws;
  p = PROTECT(coerceVector(p, REALSXP));
  m = PROTECT(coerceVector(m, REALSXP));
  draws = PROTECT(allocVector(REALSXP, n_draws));
  p_ptr = REAL(p);
  m_ptr = REAL(m);
  draws_ptr = REAL(draws);
  spline_draw(length(p), p_ptr, m_ptr, n_draws, draws_ptr);
  UNPROTECT(3);
  return draws;
}

