#include <R.h>
#include <Rinternals.h>

#ifdef __linux__
#include <malloc.h>
#endif

SEXP RAS_malloc_trim(void) {
#ifdef __linux__
    int out = malloc_trim(0);
    return ScalarInteger(out);
#else
    return ScalarInteger(NA_INTEGER);
#endif
}
