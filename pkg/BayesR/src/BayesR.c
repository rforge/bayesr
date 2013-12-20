#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for(int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}


SEXP iwls_eval(SEXP fun, SEXP response, SEXP eta)
{
  SEXP R_fcall, ans, s;
  PROTECT(R_fcall = lang3(fun, response, eta));
  s = CDR(R_fcall);
  SET_TAG(s, install("y"));
  s = CDR(R_fcall);
  SET_TAG(s, install("eta"));
  //ans = eval(R_fcall, rho);
  UNPROTECT(1);
  return R_fcall;
}


SEXP do_propose(SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id)
{
  int nProtected;

  SEXP loglik;
  PROTECT(loglik = getListElement(family, "loglik"));
  ++nProtected;
  SEXP weights;
  PROTECT(weights = getListElement(getListElement(family, "weights"), CHAR(STRING_ELT(id, 0))));
  ++nProtected;
  SEXP score;
  PROTECT(score = getListElement(getListElement(family, "score"), CHAR(STRING_ELT(id, 0))));
  ++nProtected;

  SEXP e_loglik = iwls_eval(loglik, response, eta);

  UNPROTECT(nProtected);

  return e_loglik;
}

