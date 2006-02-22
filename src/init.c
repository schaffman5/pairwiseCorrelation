#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare your C function
void Rcormat(double *x, int *nr, int *nc, double *d, int *diag, double *missingThresh);

// Define the registration table
static const R_CMethodDef CEntries[] = {
  {"Rcormat", (DL_FUNC) &Rcormat, 6}, // '6' represents number of arguments
  {NULL, NULL, 0}
};

// Define the initialization function
void R_init_pairwiseCorrelation(DllInfo *dll) {
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
