#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <math.h>

#include "multi_skew.h"
#include "skew_spline.h"
#include "skew_old.h"
#include "symmetric_Hermite.h"

void multi_skew_draw_eval(int is_draw, int n_draws, int n_grid_points, int code,
                          double mode, double *h, double mu, double omega,
                          double *z, double *ln_f)
{
  Skew_parameters skew;
  int i_draw;
  switch (code) {
    case 0: // spline draw, without v transformation
    case 1: // spline draw, with (1) v transformation
      for (i_draw = 0; i_draw < n_draws; i_draw++)
        skew_spline_draw_eval(is_draw, n_grid_points, code, mode, h, mu, omega, z+i_draw, ln_f+i_draw);
      break;

    case 2: // old skew_draw
      skew.mode = mode;
      skew.h2 = h[2] - omega; skew.h3 = h[3]; skew.h4 = h[4]; skew.h5 = h[5];
      skew.s2_prior = 1.0 / omega;
      skew.u_sign = 1.0;
      skew.is_draw = is_draw;

      /* Precomputation */
      double threshold = 0.1;
      double K_1_threshold;
      double K_2_threshold[6];
      const int max_n_reject = 100;
      Symmetric_Hermite *sh = symmetric_Hermite_alloc(100, max_n_reject); // Maximum value of K, maximum number of rejection sampling rejections
      K_1_threshold = sqrt(12.0 * threshold);
      K_2_threshold[2] = sqrt(2.0 * threshold);
      K_2_threshold[3] = exp( log(6.0 * threshold)/3.0 );
      K_2_threshold[4] = exp( log(24.0 * threshold)/24.0 );
      K_2_threshold[5] = exp( log(120.0 * threshold)/120.0 );

      for (i_draw = 0; i_draw < n_draws; i_draw++) {
        if (!is_draw)
          skew.z = z[i_draw];
        skew_draw_eval(&skew, sh, K_1_threshold, K_2_threshold);
        if (is_draw)
          z[i_draw] = skew.z;
        ln_f[i_draw] = skew.log_density;
      }
      symmetric_Hermite_free(sh);
      break;
  }
}

SEXP skew_eval_c(int n_grid_points, int code, double mode, SEXP h, double mu, double omega, SEXP z)
{
  double *h_ptr, *z_ptr, *log_f_ptr;
  SEXP log_f;
  h = PROTECT(coerceVector(h, REALSXP));
  z = PROTECT(coerceVector(z, REALSXP));
  log_f = PROTECT(allocVector(REALSXP, length(z)));
  h_ptr = REAL(h);
  z_ptr = REAL(z);
  log_f_ptr = REAL(log_f);
  multi_skew_draw_eval(FALSE, length(z), n_grid_points, code,
                       mode, h_ptr, mu, omega, z_ptr, log_f_ptr);
  UNPROTECT(3);
  return log_f;
}

SEXP skew_draw_c(int n_grid_points, int code, double mode, SEXP h, double mu, double omega, int n_draws)
{
  double *h_ptr, *draws_ptr, *log_f_ptr;
  SEXP draws, log_f;
  h = PROTECT(coerceVector(h, REALSXP));
  draws = PROTECT(allocVector(REALSXP, n_draws));
  log_f = PROTECT(allocVector(REALSXP, n_draws));
  h_ptr = REAL(h);
  draws_ptr = REAL(draws);
  log_f_ptr = REAL(log_f);
  multi_skew_draw_eval(TRUE, n_draws, n_grid_points, code,
                       mode, h_ptr, mu, omega, draws_ptr, log_f_ptr);
  UNPROTECT(3);
  return draws;
}
