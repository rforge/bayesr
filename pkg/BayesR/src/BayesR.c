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

void merr()
{
  char *m = "stopped";
  error(m);
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

SEXP do_propose(SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id, SEXP rho)
{
  int nProtected = 0;

  /* Last try accepted? */
/*  int accepted = LOGICAL(getListElement(getListElement(x, "state"), "accepted"))[0];*/
/*  int adaptive = LOGICAL(getListElement(getListElement(x, "xt"), "adaptive"))[0];*/
/*  int adaptcheck = accepted * adaptive;*/

  /* Evaluate loglik, weights and score vector. */
  SEXP eta2;
  PROTECT(eta2 = duplicate(eta));
  ++nProtected;
  int ll_ind = getListElement_index(family, "loglik");
  double pibeta = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, eta2, rho))[0];
  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, eta2, rho));
  ++nProtected;
  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, eta2, rho));
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

  /* Working observations and other stuff. */
  SEXP z, z2;
  PROTECT(z = allocVector(REALSXP, n));
  ++nProtected;
  PROTECT(z2 = allocVector(REALSXP, n));
  ++nProtected;
  SEXP state;
  PROTECT(state = duplicate(getListElement(x, "state")));
  ++nProtected;
  double *zptr = REAL(z);
  double *z2ptr = REAL(z2);
  double *etaptr = REAL(getListElement(eta2, CHAR(STRING_ELT(id, 0))));
  double *scoreptr = REAL(score);
  double *fitptr = REAL(getListElement(state, "fit"));

  /* Compute transpose of weighted design matrix, */
  /* working observations and updated predictor. */
  for(i = 0; i < n; i++) {
    for(j = 0; j < k; j++) {
      XWptr[j + k * i] = Xptr[i + n * j] * Wptr[i];
    }
    zptr[i] = etaptr[i] + scoreptr[i] / Wptr[i];
    etaptr[i] -= fitptr[i];
    z2ptr[i] = zptr[i] - etaptr[i];
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
        Pptr[i + k * j] += 1 / tau2 * Sptr[i + k * j];
      }
    }
  }

  /* Cholesky decompostion of P. */
  SEXP L;
  PROTECT(L = duplicate(P));
  ++nProtected;
  double *Lptr = REAL(L);

  for(j = 0; j < k; j++) { 	/* Zero the lower triangle. */
	  for(i = j + 1; i < k; i++) {
      Lptr[i + k * j] = 0.;
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
	F77_CALL(dpotrf)("Upper", &k, PINVLptr, &k, &info);

	for(j = 0; j < k; j++) {
	  for(i = j + 1; i < k; i++) {
		  PINVptr[i + j * k] = PINVptr[j + i * k];
    }
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
    XWptr, &k, z2ptr, &n, &zero, mu0ptr, &k);
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
  char *transa2 = "T";
  
  F77_CALL(dgemm)(transa2, transb, &k, &k1, &k, &one,
    PINVLptr, &k, g0ptr, &k, &zero, g1ptr, &k);

  double sdiag0 = 0.0;
  for(j = 0; j < k; j++) {
    g1ptr[j] += mu1ptr[j];
    sdiag0 += log(pow(Lptr[j + k * j], 2));
  }

  /* Log priors. */
  double *gptr = REAL(getListElement(state, "g"));
  double p1 = 0.0;
  double p2 = 0.0;
  double qbetaprop = 0.0;

  double b = 0.0;
  double tsum1 = 0.0;
  double tsum2 = 0.0;
  double tsum3 = 0.0;
  double tsum4 = 0.0;
  for(i = 0; i < k; i++) {
    tsum1 = 0.0;
    tsum2 = 0.0;
    tsum3 = 0.0;
    tsum4 = 0.0;
    for(j = 0; j < k; j++) {
      if(fixed < 1) {
        tsum1 += g1ptr[j] * Sptr[j + i * k];
        tsum2 += gptr[j] * Sptr[j + i * k];
        tsum4 += g1ptr[j] * Pptr[j + i * k];
      }
      tsum3 += (g1ptr[j] - mu1ptr[j]) * Pptr[j + i * k];
    }
    if(fixed < 1) {
      p1 += tsum2 * gptr[i];
      p2 += tsum1 * g1ptr[i];
      b += tsum4 * g1ptr[i];
    }
    qbetaprop += tsum3 * (g1ptr[i] - mu1ptr[i]);
  }

  if(fixed < 1) {
    p1 = -0.5 / tau2 * p1;
    p2 = -0.5 / tau2 * p2;
  }

  qbetaprop = 0.5 * sdiag0 - 0.5 * qbetaprop;

  /* Obtain new fitted values and update predictor. */
  SEXP fit1;
  PROTECT(fit1 = allocVector(REALSXP, n));
  ++nProtected;
  double *fit1ptr = REAL(fit1);

  for(i = 0; i < n; i++) {
    fit1ptr[i] = 0.0;
    for(j = 0; j < k; j++) {
      fit1ptr[i] += Xptr[i + n * j] * g1ptr[j];
    }
    etaptr[i] = etaptr[i] + fit1ptr[i];
  }

  double pibetaprop = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, eta2, rho))[0];

  /* Weights, score and working observations. */
  SEXP weights2, score2;
  PROTECT(weights2 = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, eta2, rho));
  ++nProtected;
  PROTECT(score2 = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, eta2, rho));
  ++nProtected;
  double *W2ptr = REAL(weights2);
  double *score2ptr = REAL(score2);

  /* Compute transpose of weighted design matrix */
  for(i = 0; i < n; i++) {
    for(j = 0; j < k; j++) {
      XWptr[j + k * i] = Xptr[i + n * j] * W2ptr[i];
    }
    zptr[i] = etaptr[i] + score2ptr[i] / W2ptr[i];
    etaptr[i] -= fit1ptr[i];
    z2ptr[i] = zptr[i] - etaptr[i];
  }

  /* Compute X'WX. */
  F77_CALL(dgemm)(transa, transb, &k, &k, &n, &one,
    XWptr, &k, Xptr, &n, &zero, Pptr, &k);

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(i = 0; i < k; i++) {
      for(j = 0; j < k; j++) {
        Pptr[i + k * j] += 1 / tau2 * Sptr[i + k * j];
      }
    }
  }

  /* Cholesky decompostion of P. */
  L = duplicate(P);
  Lptr = REAL(L);
  for(j = 0; j < k; j++) { 	/* Zero the lower triangle. */
	  for(i = j + 1; i < k; i++) {
      Lptr[i + k * j] = 0.;
    }
  }

	F77_CALL(dpotrf)("Upper", &k, Lptr, &k, &info);

  /* Compute the inverse precision matrix. */
  PINV = duplicate(L);
  PINVptr = REAL(PINV);
  F77_CALL(dpotri)("Upper", &k, PINVptr, &k, &info);

	for(j = 0; j < k; j++) {
	  for(i = j + 1; i < k; i++) {
		  PINVptr[i + j * k] = PINVptr[j + i * k];
    }
  }

  /* Compute mu. */
  F77_CALL(dgemm)(transa, transb, &k, &k1, &n, &one,
    XWptr, &k, z2ptr, &n, &zero, mu0ptr, &k);
  F77_CALL(dgemm)(transa, transb, &k, &k1, &k, &one,
    PINVptr, &k, mu0ptr, &k, &zero, mu1ptr, &k);

  double qbeta = 0.0;
  sdiag0 = 0.0;

  for(i = 0; i < k; i++) {
    tsum3 = 0.0;
    for(j = 0; j < k; j++) {
      tsum3 += (gptr[j] - mu1ptr[j]) * Pptr[j + i * k];
    }
    qbeta += tsum3 * (gptr[i] - mu1ptr[i]);
    sdiag0 += log(pow(Lptr[i + k * i], 2));
  }

  qbeta = 0.5 * sdiag0 - 0.5 * qbeta;

  SEXP tau3;
  PROTECT(tau3 = allocVector(REALSXP, 1));
  ++nProtected;
  if(fixed < 1) {
    b += REAL(getListElement(x, "b"))[0];
    GetRNGstate();
    REAL(tau3)[0] = 1 / rgamma(REAL(getListElement(x, "rank"))[0], 1 / b);
    PutRNGstate();
  } else {
    REAL(tau3)[0] = tau2;
  }

  SEXP alpha;
  PROTECT(alpha = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(alpha)[0] = (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1);

  /* Stuff everything together. */
  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 4));
  ++nProtected;

  SET_VECTOR_ELT(rval, 0, g1);
  SET_VECTOR_ELT(rval, 1, fit1);
  SET_VECTOR_ELT(rval, 2, tau3);
  SET_VECTOR_ELT(rval, 3, alpha);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 4));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("g"));
  SET_STRING_ELT(nrval, 1, mkChar("fit"));
  SET_STRING_ELT(nrval, 2, mkChar("tau2"));
  SET_STRING_ELT(nrval, 3, mkChar("alpha"));
        
  setAttrib(rval, R_NamesSymbol, nrval); 

  UNPROTECT(nProtected);

  return rval;
}

