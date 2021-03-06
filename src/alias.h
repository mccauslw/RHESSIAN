#ifndef ALIAS
#define ALIAS

#include <R.h>
#include <Rinternals.h>

// For internal C use
void draw_discrete(int n, double *p, int n_draws, int *draws);
void alias_tables(int n, double *p, int *Alias, double *Prob);
void draw_discrete_from_alias_tables(int n, int *Alias, double *Prob, int n_draws, int *draws);

// For use in C++ wrapper used for testing
SEXP alias_c(SEXP p, int n_draws);

#endif /* ALIAS */
