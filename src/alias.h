#ifndef ALIAS
#define ALIAS

// For internal C use
void draw_discrete(int n, double *p, int n_draws, int *draws);

// For use in C++ wrapper used for testing
SEXP alias_c(SEXP p, int n_draws);

#endif /* ALIAS */
