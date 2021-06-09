#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#define PHY_API_IMPLEMENTATION
#include <phy.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


SEXP parsimony_fitch_mpr(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP parsimony_fitch_mpr_count(SEXP, SEXP, SEXP, SEXP);
SEXP parsimony_fitch_mpr_sample(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP parsimony_sankoff_mpr_downpass(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP parsimony_sankoff_mpr_uppass(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP parsimony_sankoff_mpr_count(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP parsimony_sankoff_mpr_sample(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    CALLDEF(parsimony_fitch_mpr, 5),
    CALLDEF(parsimony_fitch_mpr_count, 4),
    CALLDEF(parsimony_fitch_mpr_sample, 5),
    CALLDEF(parsimony_sankoff_mpr_downpass, 5),
    CALLDEF(parsimony_sankoff_mpr_uppass, 6),
    CALLDEF(parsimony_sankoff_mpr_count, 6),
    CALLDEF(parsimony_sankoff_mpr_sample, 8),
    {NULL, NULL, 0}
};


void attribute_visible R_init_parsimony(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
