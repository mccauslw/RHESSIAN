#include <R.h>
#include <Rinternals.h>
#include "alias.h"
#include "RNG.h"
#include "Phi.h"

// Vose's Alias Method to draw a discrete random variable
#define max_n 100
void draw_discrete(int n, double *p, int n_draws, int *draws)
{
  int i, l, g;
  int i_draw;
  double sum = 0.0;
  double P[max_n]; // Normalized probabilities times n
  rng_init_rand(5489UL);

  // 1. Create arrays alias and prob
  int Alias[max_n];
  double Prob[max_n];

  // 2. Create two worklists, Small and Large
  int Large[max_n]; int n_Large = 0;
  int Small[max_n]; int n_Small = 0;

  // 3. Multiply each probability by n
  for (i=0; i<n; i++)
    sum += p[i];
  for (i=0; i<n; i++)
    P[i] = p[i] * n / sum;

  // 4. For each scaled probability P[i]:
  //    (a) if P[i] < 1, add i to Small;
  //    (b) otherwise add i to Large
  for (i=0; i<n; i++)
    if (P[i] < 1) Small[n_Small++] = i; else Large[n_Large++] = i;

  // 5. While Small and Large are not empty:
  for (; n_Large > 0 && n_Small > 0; ) {
    // (a) Remove first element from Small, call it l
    l = Small[--n_Small];
    // (b) Remove first element from Large, call it g
    g = Large[--n_Large];
    // (c) and (d) Set prob[l] = p[l], alias[l] = g
    Prob[l] = P[l];
    Alias[l] = g;
    // (e) Set p_g = P[g] + P[l] - 1
    P[g] = P[g] + P[l] - 1.0;
    // (f) and (g) If p_g < 1, add g to small; otherwise add g to Large
    if (P[g] < 1) Small[n_Small++] = g; else Large[n_Large++] = g;
  }

  // 6. While Large is not empty: remove first element from Large, call it g, set P[g] = 1
  for (; n_Large > 0; ) {
    P[Large[--n_Large]] = 1.0;
    Alias[Large[n_Large]] = Large[n_Large];
  }

  // 7. While Small is not empty: remove first element from Small, call it l, set P[l] = 1
  for (; n_Small > 0; ) {
    P[Small[--n_Small]] = 1.0;
    Alias[Small[n_Small]] = Small[n_Small];
  }

  // Generation
  for (i_draw=0; i_draw<n_draws; i_draw++) {
    i = floor(n * rng_rand());
    draws[i_draw] = (rng_rand() < Prob[i]) ? i : Alias[i];
  }
}

SEXP alias_c(SEXP p, int n_draws) {
  double *p_ptr;
  int *draws_ptr;
  SEXP draws;
  p = PROTECT(coerceVector(p, REALSXP));
  draws = PROTECT(allocVector(INTSXP, n_draws));
  p_ptr = REAL(p);
  draws_ptr = INTEGER(draws);
  draw_discrete(length(p), p_ptr, n_draws, draws_ptr);
  UNPROTECT(2);
  return draws;
}
