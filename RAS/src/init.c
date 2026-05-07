#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP RAS_malloc_trim(void);

static const R_CallMethodDef CallEntries[] = {
    {"RAS_malloc_trim", (DL_FUNC) &RAS_malloc_trim, 0},
    {NULL, NULL, 0}
};

void R_init_RAS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
