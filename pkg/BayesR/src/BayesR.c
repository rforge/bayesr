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


SEXP iwls_eval(SEXP fun, SEXP response, SEXP eta, SEXP rho)
{
  SEXP R_fcall, args, rval;

  PROTECT(R_fcall = lang3(fun, response, eta));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


SEXP do_propose(SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id, SEXP rho)
{
  int nProtected;

  /* Obtain loglik, weights and score function. */
  SEXP loglik, weights, score, e_weights, e_score;
  PROTECT(loglik = getListElement(family, "loglik"));
  ++nProtected;
  PROTECT(weights = getListElement(getListElement(family, "weights"), CHAR(STRING_ELT(id, 0))));
  ++nProtected;
  PROTECT(score = getListElement(getListElement(family, "score"), CHAR(STRING_ELT(id, 0))));
  ++nProtected;

  /* Evaluate loglik, weights and score vector. */
  double pibeta = REAL(iwls_eval(loglik, response, eta, rho))[0];
  PROTECT(e_weights = iwls_eval(weights, response, eta, rho));
  ++nProtected;
  PROTECT(e_score = iwls_eval(score, response, eta, rho));
  ++nProtected;

  /* Extract design matrix X and penalty matrix S.*/
  SEXP X, S;
  PROTECT(X = getListElement(x, "X"));
  ++nProtected;
  PROTECT(S = VECTOR_ELT(getListElement(x, "S"), 0));
  ++nProtected;

  int i, j, n = nrows(X), k = ncols(S);
  double *Xptr = REAL(X);
  double *Sptr = REAL(S);
  double *Wptr = REAL(e_weights);
  double tau2 = REAL(getListElement(getListElement(x, "state"), "tau2"))[0];

  /* Create weighted matrix */
  SEXP XW, dim;
  PROTECT(XW = allocVector(REALSXP, k * n));
  ++nProtected;
  PROTECT(dim = allocVector(INTSXP, 2));
  ++nProtected;
  INTEGER(dim)[0] = k; 
  INTEGER(dim)[1] = n;   
  setAttrib(XW, R_DimSymbol, dim);
  double *XWptr = REAL(XW);

  /* Compute transpose of weighted design matrix */
  for(j = 0; j < k; ++i) {
    for(i = 0; i < n; ++j) {
      XWptr[j * n + i] = 1; // Xptr[j * n + i] * Wptr[i];
    }
  }

Rprintf("ok\n");

  UNPROTECT(nProtected);

  return XW;
}

