#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


SEXP test_fc(SEXP strings, SEXP refstring);
//SEXP additiveFunctionValues(SEXP reachIDs, SEXP binaryIDs, SEXP proportionalInfluences, SEXP netIDs);
//SEXP DoWritedbfSSN(SEXP file, SEXP df, SEXP pr, SEXP sc, SEXP DataTypes);

static const R_CallMethodDef CallEntries[] = {
	CALLDEF(test_fc, 2),
	//CALLDEF(additiveFunctionValues, 4),
	//CALLDEF(DoWritedbfSSN, 5),
    {NULL, NULL, 0}
};
void R_init_SSN2test(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
