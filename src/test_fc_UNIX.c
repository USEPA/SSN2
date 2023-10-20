/*
  test_fc.c:
    Compare a vector of strings to a reference string. For each string
    in the vector, return the (0-based) index of the first character
    that differs from the reference string, which in R equates to the
    (1-based) index of the last matching character. If the strings are
    identical up to the size of the shorter string, return -1 times the
    size. 

    The integer vector returned to R will contain positive integers
    corresponding to the number of character

    example:
      test_fc(c("1", "111", "101", "011"), "11")  -->  c(-1, -2, 1, 0)

  Jim Regetz
  NCEAS
*/

#include <R.h>
#include <Rinternals.h>
SEXP test_fc(SEXP strings, SEXP refstring) {
  SEXP ans;
  int i, j, len;
  int nchar1, nchar2, nchar;
  const char *ref;

  // lightweight error checking
  if(!isString(strings) || !isString(refstring)) {
    error("invalid arguments");
  }
  if(length(refstring)!=1) {
    error("second argument must be a length 1 character vector");
  }

  // extract reference string
  ref = CHAR(STRING_ELT(refstring, 0));
  nchar2 = strlen(ref);

  len = length(strings);
  PROTECT(ans = allocVector(INTSXP, len));
  for(i = 0; i < len; i++) {

    nchar1 = strlen(CHAR(STRING_ELT(strings, i)));
    nchar = nchar1<nchar2 ? nchar1 : nchar2;
    INTEGER(ans)[i] = -1*nchar;
    // record index of first differing character, if any
    for (j=0; j<nchar; j++) {
      if (CHAR(STRING_ELT(strings,i))[j]!=ref[j]) {
        INTEGER(ans)[i] = j;
        break;
      }
    }
  }
  UNPROTECT(1);
  return ans;
}
