#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include <R_ext/Applic.h> /* for dgemm */
#include <R_ext/Complex.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Linpack.h>


/* (1) Helper functions. */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt, names;
  PROTECT(elmt = R_NilValue);
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  for(int i = 0; i < length(list); i++) {
	  if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	  }
  }

  UNPROTECT(2);

  return elmt;
}


int getListElement_index(SEXP list, const char *str)
{
  SEXP names;
  PROTECT(names = getAttrib(list, R_NamesSymbol));

  int j = -1;

  for(int i = 0; i < length(list); i++) {
	  if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    j = i;
	    break;
	  }
  }

  UNPROTECT(1);

  return j;
}


void pvec(SEXP vec)
{
  int i, n = length(vec);
  double *vecptr = REAL(vec);
  for(i = 0; i < n; i++, vecptr++)
    Rprintf(" %g", *vecptr);
  Rprintf("\n");
}
    
void pmat(SEXP mat)
{
  int i,j,n = nrows(mat), k = ncols(mat);
  Rprintf("   ");
  for(j = 0; j < k; ++j)
    Rprintf("[%d] ", j);
  Rprintf("\n");
  for(i = 0; i < n; ++i) {
    Rprintf("[%d]", i);
    for(j = 0; j < k; ++j)
      Rprintf(" %g", REAL(mat)[i + j * n]);
    Rprintf("\n");
  }
  Rprintf("\n");
}


/* Matrix product. */
SEXP matprod(SEXP x, SEXP y)
{
  int nrx, ncx, nry, ncy;
  double *ansptr, *xptr, *yptr;

  nrx = nrows(x);
  ncx = ncols(x);
  nry = nrows(y);
  ncy = ncols(y);

  //PROTECT(x = coerceVector(x, REALSXP));
  //PROTECT(y = coerceVector(y, REALSXP));
  xptr = REAL(x);  yptr = REAL(y);

  SEXP ans;
  PROTECT(ans = allocMatrix(REALSXP, nrx, ncy));
  ansptr = REAL(ans);

  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
    xptr, &nrx, yptr, &nry, &zero, ansptr, &nrx);

  UNPROTECT(1);

  return(ans);
}


/* (2) Main IWLS functions. */
SEXP iwls_eval(SEXP fun, SEXP response, SEXP eta, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang3(fun, response, eta));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}

/*https://gist.github.com/Sharpie/323498*/
SEXP do_propose(SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id, SEXP rho)
{
  int nProtected;

  /* Evaluate loglik, weights and score vector. */
  double pibeta = REAL(iwls_eval(getListElement(family, "loglik"), response, eta, rho))[0];
  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, eta, rho));
  ++nProtected;
  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, eta, rho));
  ++nProtected;

  /* Extract design matrix X and penalty matrix S.*/
  int X_ind = getListElement_index(x, "X");
  int S_ind = getListElement_index(x, "S");
  double *Xptr = REAL(VECTOR_ELT(x, X_ind));
  double *Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), 0));
  double *Wptr = REAL(weights);

  int i, j, n = nrows(VECTOR_ELT(x, X_ind)), k = ncols(VECTOR_ELT(x, X_ind));
  double tau2 = REAL(getListElement(getListElement(x, "state"), "tau2"))[0];

  /* Create weighted matrix */
  SEXP XW;
  PROTECT(XW = allocMatrix(REALSXP, k, n));
  ++nProtected;
  double *XWptr = REAL(XW);

  /* Compute transpose of weighted design matrix */
  for(i = 0; i < n; i++) {
    for(j = 0; j < k; j++) {
      XWptr[j + k * i] = Xptr[i + n * j] * Wptr[i];
    }
  }

  /* Compute X'WX. */
  SEXP P;
  PROTECT(P = allocMatrix(REALSXP, k, k));
  ++nProtected;
  double *Pptr = REAL(P);
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, &k, &k, &n, &one,
    XWptr, &k, Xptr, &n, &zero, Pptr, &k);

  /* Add penalty matrix and variance parameter. */
  int fixed = LOGICAL(getListElement(x, "fixed"))[0];
  if(fixed < 1) {
    for(i = 0; i < k; i++) {
      for(j = 0; j < k; j++) {
        Pptr[i + k * j] += tau2 * Sptr[i + k * j];
      }
    }
  }

  /* Cholesky decompostion of P. */
  SEXP L;
  PROTECT(L = allocMatrix(REALSXP, k, k));
  ++nProtected;
  double *Lptr = REAL(L);

  for(i = 0; i < k; i++) {
	  for(j = 0; j < k; j++) {
      if(i > j)
        Lptr[i + k * j] = 0.0;
      else
        Lptr[i + k * j] = Pptr[i + k * j];
    }
  }

	int info;
	F77_CALL(dpotrf)("Upper", &k, Lptr, &k, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV;
  PROTECT(PINV = duplicate(L));
  ++nProtected;
  double *PINVptr = REAL(PINV);

  F77_CALL(dpotri)("Upper", &k, PINVptr, &k, &info);

  SEXP PINVL;
  PROTECT(PINVL = duplicate(PINV));
  ++nProtected;
  double *PINVLptr = REAL(PINVL);

	for(j = 0; j < k; j++) {
	  for(i = j + 1; i < k; i++) {
		  PINVptr[i + j * k] = PINVptr[j + i * k];
    }
  }

  /* Compute the working observations. */
  SEXP z;
  PROTECT(z = allocVector(REALSXP, n));
  ++nProtected;
  SEXP state;
  PROTECT(state = duplicate(getListElement(x, "state")));
  ++nProtected;
  double *zptr = REAL(z);
  double *etaptr = REAL(getListElement(eta, CHAR(STRING_ELT(id, 0))));
  double *scoreptr = REAL(score);
  double *fitptr = REAL(getListElement(state, "fit"));

  for(i = 0; i < n; i++) {
    zptr[i] = scoreptr[i] / Wptr[i] + fitptr[i];
  }

  /* Compute mu. */
  int k1 = 1;
  SEXP mu0, mu1;
  PROTECT(mu0 = allocVector(REALSXP, k));
  ++nProtected;
  PROTECT(mu1 = allocVector(REALSXP, k));
  ++nProtected;
  double *mu0ptr = REAL(mu0);
  double *mu1ptr = REAL(mu1);
  F77_CALL(dgemm)(transa, transb, &k, &k1, &n, &one,
    XWptr, &k, zptr, &n, &zero, mu0ptr, &k);
  F77_CALL(dgemm)(transa, transb, &k, &k1, &k, &one,
    PINVptr, &k, mu0ptr, &k, &zero, mu1ptr, &k);

  /* Sample. */
  SEXP g0;
  PROTECT(g0 = allocVector(REALSXP, k));
  ++nProtected;
  double *g0ptr = REAL(g0);

  GetRNGstate();
  for(j = 0; j < k; j++) {
    g0ptr[j] = rnorm(0, 1);
  }
  PutRNGstate();

  SEXP g1;
  PROTECT(g1 = allocVector(REALSXP, k));
  ++nProtected;
  double *g1ptr = REAL(g1);

  F77_CALL(dgemm)(transa, transb, &k, &k1, &k, &one,
    PINVLptr, &k, g0ptr, &k, &zero, g1ptr, &k);

  for(j = 0; j < k; j++) {
    g1ptr[j] += mu1ptr[j];
  }

  /* Log priors. */
  double *gptr = REAL(getListElement(state, "g"));
  double p1 = 0.0;
  double p2 = 0.0;
  if(fixed < 1) {
    double tsum = 0.0;
    double tsum2 = 0.0;
    for(i = 0; i < k; i++) {
      tsum = 0.0;
      tsum2 = 0.0;
      for(j = 0; j < k; j++) {
        tsum += g1ptr[j] * Sptr[j + i * k];
        tsum2 += gptr[j] * Sptr[j + i * k];
      }
      p1 += tsum2 * gptr[i];
      p2 += tsum * g1ptr[i];
    }
  }

  // qbeta <- 0.5*sum(log((diag(cholprprop)^2)))-0.5*t(coeff-muprop)%*%prprop%*%(coeff-muprop)

  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 4));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, P);
  SET_VECTOR_ELT(rval, 1, L);
  SET_VECTOR_ELT(rval, 2, g1);
  SET_VECTOR_ELT(rval, 3, mu1);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 4));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("P"));
  SET_STRING_ELT(nrval, 1, mkChar("L"));
  SET_STRING_ELT(nrval, 2, mkChar("g"));
  SET_STRING_ELT(nrval, 3, mkChar("mu"));

  setAttrib(rval, R_NamesSymbol, nrval);        

  UNPROTECT(nProtected);

  return rval;
}

