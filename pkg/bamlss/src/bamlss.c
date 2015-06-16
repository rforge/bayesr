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


/* Map linear predictor to parameter scale */
SEXP map2par(SEXP fun, SEXP eta, SEXP rho)
{
  SEXP R_fcall, rval;

  PROTECT(R_fcall = lang2(fun, eta));
  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);

  return rval;
}


SEXP do_propose(SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id, SEXP rho)
{
  int nProtected = 0;
  int fixed = LOGICAL(getListElement(x, "fixed"))[0];
  int fxsp = LOGICAL(getListElement(x, "fxsp"))[0];

  /* Last try accepted? */
/*  int accepted = LOGICAL(getListElement(getListElement(x, "state"), "accepted"))[0];*/
/*  int adaptive = LOGICAL(getListElement(getListElement(x, "xt"), "adaptive"))[0];*/
/*  int adaptcheck = accepted * adaptive;*/

  /* Evaluate loglik, weights and score vector. */
  SEXP eta2, peta;
  PROTECT(eta2 = duplicate(eta));
  ++nProtected;
  PROTECT(peta = map2par(getListElement(family, "map2par"), eta2, rho));
  ++nProtected;
  int ll_ind = getListElement_index(family, "loglik");
  double pibeta = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, rho))[0];
  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho));
  ++nProtected;
  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho));
  ++nProtected;

  /* Extract design matrix X and penalty matrix S.*/
  int X_ind = getListElement_index(x, "X");
  double *Xptr = REAL(VECTOR_ELT(x, X_ind));
  double *Wptr = REAL(weights);
  int S_ind = 0;
  double * Sptr;
  if(fixed < 1) {
    S_ind = getListElement_index(x, "S");
    Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), 0));
  }

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
    sdiag0 += log(pow(Lptr[j + k * j], 2.0));
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

  peta = map2par(getListElement(family, "map2par"), eta2, rho);
  double pibetaprop = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, rho))[0];

  /* Weights, score and working observations. */
  SEXP weights2, score2;
  PROTECT(weights2 = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho));
  ++nProtected;
  PROTECT(score2 = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho));
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
    sdiag0 += log(pow(Lptr[i + k * i], 2.0));
  }

  qbeta = 0.5 * sdiag0 - 0.5 * qbeta;

  SEXP tau3;
  PROTECT(tau3 = allocVector(REALSXP, 1));
  ++nProtected;
  if(fixed < 1 && fxsp < 1) {
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


/* logPost evaluation */
double lp_eval(SEXP fun, SEXP g, SEXP x,
  SEXP family, SEXP response, SEXP eta,
  SEXP id, SEXP rho)
{
  SEXP R_fcall, t, rval;

  PROTECT(t = R_fcall = allocList(7));
  SET_TYPEOF(R_fcall, LANGSXP);

  SETCAR(R_fcall, fun);
  t = CDR(t); SETCAR(t, g);
  t = CDR(t); SETCAR(t, x);
  t = CDR(t); SETCAR(t, family);
  t = CDR(t); SETCAR(t, response);
  t = CDR(t); SETCAR(t, eta);
  t = CDR(t); SETCAR(t, id);

  PROTECT(rval = eval(R_fcall, rho));

  UNPROTECT(2);
  return REAL(rval)[0];
}


/* Univariate slice sampling */
SEXP uni_slice(SEXP g, SEXP x, SEXP family, SEXP response, SEXP eta, SEXP id, SEXP j,
  SEXP W, SEXP M, SEXP LOWER, SEXP UPPER, SEXP logPost, SEXP rho)
{
  int nProtected = 0;
  int jj = INTEGER(j)[0] - 1;

  int m = INTEGER(M)[0] + 1;
  double w = REAL(W)[0];
  double lower = REAL(LOWER)[0];
  double upper = REAL(UPPER)[0];

  SEXP gL, gR;
  PROTECT(gL = duplicate(g));
  ++nProtected;
  PROTECT(gR = duplicate(g));
  ++nProtected;

  double *gLptr = REAL(gL);
  double *gRptr = REAL(gR);
  double *gptr = REAL(g);

  double x0 = gptr[jj];
  double gx0 = lp_eval(logPost, g, x, family, response, eta, id, rho);

  GetRNGstate();
  double logy = gx0 - rexp(1); 
  double u = runif(0.0, w);
  PutRNGstate();

  gLptr[jj] = x0 - u;
  gRptr[jj] = x0 + (w - u);

  if(m > 1) {
    GetRNGstate();
    int J = floor(runif(0.0, m));
    PutRNGstate();
    int K = (m - 1) - J;
    while(J > 0) {
      if(gLptr[jj] <= lower)
        break;
      if(lp_eval(logPost, gL, x, family, response, eta, id, rho) <= logy)
        break;
      gLptr[jj] = gLptr[jj] - w;
      J = J - 1;
    }
    while(K > 0) {
      if(gRptr[jj] >= upper)
        break;
      if(lp_eval(logPost, gR, x, family, response, eta, id, rho) <= logy)
        break;
      gRptr[jj] = gRptr[jj] + w;
      K = K - 1;
    }
  }

  if(gLptr[jj] < lower) {
    gLptr[jj] = lower;
  }
  if(gRptr[jj] > upper) {
    gRptr[jj] = upper;
  }

  int run = 1;
  
  while(run > 0) {
    gptr[jj] = runif(gLptr[jj], gRptr[jj]);

    double gx1 = lp_eval(logPost, g, x, family, response, eta, id, rho);

    if(gx1 >= logy)
      run = 0;

    if(gptr[jj] > x0) {
      gRptr[jj] = gptr[jj];
    } else {
      gLptr[jj] = gptr[jj];
    }
  }

  UNPROTECT(nProtected);
  return g;
}


/* Compute the centroid of a polygon. */
SEXP cpos(SEXP p, SEXP K, SEXP pos)
{
  int i, n, k;
  n = INTEGER(K)[0];
  k = n + 1;
    
  double tmp, *pptr, asum, xsum, ysum;
    
  pptr = REAL(p);
    
  asum = 0.0;
  xsum = 0.0;
  ysum = 0.0;
    
  for(i = 0; i < n; i++) {
    tmp = pptr[i] * pptr[i + k + 1] - pptr[i + 1] * pptr[i + k];
    asum += tmp;
    xsum += (pptr[i] + pptr[i + 1]) * tmp;
    ysum += (pptr[i + k] + pptr[i + k + 1]) * tmp;
  }
        
  tmp = 1 / (3 * asum);
  REAL(pos)[0] = tmp * xsum;
  REAL(pos)[1] = tmp * ysum;
    
  return pos;
}


/* Compute reduced weights and residuals. */
void xbin_fun(SEXP ind, SEXP weights, SEXP e, SEXP xweights, SEXP xrres, SEXP order)
{
  int i;
  int j = 0;
  int n = length(ind);
  int k = 0;

  double *weightsptr = REAL(weights);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(xweights);
  double *xrresptr = REAL(xrres);
  int *indptr = INTEGER(ind);
  int *orderptr = INTEGER(order);

  xweightsptr[0] = 0;
  xrresptr[0] = 0;

  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      ++j;
      xweightsptr[j] = 0;
      xrresptr[j] = 0;
    }
    k = orderptr[i] - 1;
    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }
}


/* Efficient IWLS sampling. */
SEXP gmcmc_iwls(SEXP family, SEXP theta, SEXP id,
  SEXP eta, SEXP response, SEXP x, SEXP z, SEXP e, SEXP rho)
{
  int i, j, k, nProtected = 0;
  int n = INTEGER(getListElement(x, "nobs"))[0];
  int fixed = LOGICAL(getListElement(x, "fixed"))[0];
  int fxsp = LOGICAL(getListElement(x, "fxsp"))[0];

  SEXP theta2;
  PROTECT(theta2 = duplicate(getListElement(getListElement(theta,
    CHAR(STRING_ELT(id, 0))), CHAR(STRING_ELT(id, 1)))));
  ++nProtected;

  int S_ind = getListElement_index(x, "S");
  int ntau2;
  if(fixed > 0) {
    ntau2 = 0;
  } else {
    ntau2 = length(VECTOR_ELT(x, S_ind));
  }
  int nc = length(theta2);
  if(fixed < 1) {
    nc -= ntau2;
  }

  SEXP gamma0, gamma1, tau2;
  PROTECT(gamma0 = allocVector(REALSXP, nc));
  ++nProtected;
  PROTECT(gamma1 = allocVector(REALSXP, nc));
  ++nProtected;
  PROTECT(tau2 = allocVector(REALSXP, ntau2));
  ++nProtected;

  /* Evaluate loglik, weights and score vector. */
  SEXP peta;
  PROTECT(peta = map2par(getListElement(family, "map2par"), eta, rho));
  ++nProtected;
  int ll_ind = getListElement_index(family, "loglik");
  double pibeta = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, rho))[0];
  SEXP weights;
  PROTECT(weights = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho));
  double *weightsptr = REAL(weights);
  ++nProtected;
  SEXP score;
  PROTECT(score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho));
  double *scoreptr = REAL(score);
  ++nProtected;
  SEXP eta2;
  PROTECT(eta2 = duplicate(eta));
  ++nProtected;

  /* Create weighted matrix */
  int X_ind = getListElement_index(x, "X");
  int nr = nrows(VECTOR_ELT(x, X_ind));

  /* More pointers needed. */
  double *thetaptr = REAL(theta2);
  double *gamma0ptr = REAL(gamma0);
  double *gamma1ptr = REAL(gamma1);
  double *tau2ptr = REAL(tau2);
  double *etaptr = REAL(getListElement(eta2, CHAR(STRING_ELT(id, 0))));
  double *zptr = REAL(z);
  double *eptr = REAL(e);
  double *xweightsptr = REAL(getListElement(x, "weights"));
  double *xrresptr = REAL(getListElement(x, "rres"));
  double *XWptr = REAL(getListElement(x, "XW"));
  double *XWXptr = REAL(getListElement(x, "XWX"));
  double *Xptr = REAL(VECTOR_ELT(x, X_ind));
  double *Sptr;
  int *idptr = INTEGER(getListElement(x, "xbin.ind"));
  int *indptr = INTEGER(getListElement(x, "xbin.sind"));
  int *orderptr = INTEGER(getListElement(x, "xbin.order"));

  /* Handling fitted.values. */
  SEXP fitname;
  PROTECT(fitname = allocVector(STRSXP, 1));
  ++nProtected;
  SET_STRING_ELT(fitname, 0, mkChar("fitted.values"));
  double *fitptr = REAL(getAttrib(theta2, fitname));
  double *fitrptr = REAL(getListElement(x, "fit.reduced"));

  /* Start. */
  xweightsptr[0] = 0;
  xrresptr[0] = 0;

  j = 0;
  int jj;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0;
      xrresptr[j] = 0;
    }
    k = orderptr[i] - 1;

    zptr[k] = etaptr[k] + scoreptr[k] / weightsptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];

    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  char *transa = "N", *transb = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      tau2ptr[jj] = thetaptr[nc + jj];
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L;
  PROTECT(L = duplicate(getListElement(x, "XWX")));
  ++nProtected;
  double *Lptr = REAL(L);

  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
	  for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

	int info;
	F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV;
  PROTECT(PINV = duplicate(L));
  ++nProtected;
  double *PINVptr = REAL(PINV);

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);

  SEXP PINVL;
  PROTECT(PINVL = duplicate(PINV));
  ++nProtected;
  double *PINVLptr = REAL(PINVL);
	F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

	for(j = 0; j < nc; j++) {
	  for(i = j + 1; i < nc; i++) {
		  PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  int k1 = 1;
  SEXP mu0, mu1;
  PROTECT(mu0 = allocVector(REALSXP, nc));
  ++nProtected;
  PROTECT(mu1 = allocVector(REALSXP, nc));
  ++nProtected;
  double *mu0ptr = REAL(mu0);
  double *mu1ptr = REAL(mu1);
  char *transa2 = "T";
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);
  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  /* Sample. */
  GetRNGstate();
  for(j = 0; j < nc; j++) {
    gamma0ptr[j] = rnorm(0, 1);
  }
  PutRNGstate();
  
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nc, &one,
    PINVLptr, &nc, gamma0ptr, &nc, &zero, gamma1ptr, &nc);

  double sdiag0 = 0.0;
  for(j = 0; j < nc; j++) {
    gamma1ptr[j] += mu1ptr[j];
    sdiag0 += log(pow(Lptr[j + nc * j], 2.0));
  }

  /* Log priors. */
  double p2 = 0.0;
  double qbetaprop = 0.0;
  double tsum1 = 0.0;
  double tsum2 = 0.0;
  double gSg = 0.0;
  double a = 0.0;
  double b = 0.0;
  double *rankptr;
  if(fixed < 1) {
    a = REAL(getListElement(x, "a"))[0];
    b = REAL(getListElement(x, "b"))[0];
    rankptr = REAL(getListElement(x, "rank"));
    for(jj = 0; jj < ntau2; jj++) {
      Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      gSg = 0.0;
      for(i = 0; i < nc; i++) {
        tsum1 = 0.0;
        tsum2 = 0.0;
        for(j = 0; j < nc; j++) {
          tsum1 += gamma1ptr[j] * Sptr[j + i * nc];
          tsum2 += (gamma1ptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
        }
        gSg += tsum1 * gamma1ptr[i];
        qbetaprop += tsum2 * (gamma1ptr[i] - mu1ptr[i]);
      }
      p2 += -log(tau2ptr[jj]) * rankptr[jj] / 2 - 0.5 / tau2ptr[jj] * gSg +
        log(pow(b, a)) - exp(lgamma(a)) + (-a - 1) * log(tau2ptr[jj]) - b / tau2ptr[jj];
    }
  } else {
    for(i = 0; i < nc; i++) {
      tsum1 = 0.0;
      p2 += dnorm(gamma1ptr[i], 0, 1000, 1);
      for(j = 0; j < nc; j++) {
        tsum1 += (gamma1ptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
      }
      qbetaprop += tsum1 * (gamma1ptr[i] - mu1ptr[i]);
    }
  }
  qbetaprop = 0.5 * sdiag0 - 0.5 * qbetaprop;

  /* Compute edf. */
/*  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,*/
/*    XWXptr, &nc, PINVptr, &nc, &zero, mu1ptr, &nc);*/

  /* Part 2. */
  /* Obtain new fitted values and update predictor. */
  F77_CALL(dgemm)(transa, transb, &nr, &k1, &nc, &one,
    Xptr, &nr, gamma1ptr, &nr, &zero, fitrptr, &nr);

/*  for(i = 0; i < nr; i++) {*/
/*    fitrptr[i] = 0.0;*/
/*    for(j = 0; j < nc; j++) {*/
/*      fitrptr[i] += Xptr[i + nr * j] * gamma1ptr[j];*/
/*    }*/
/*  }*/

  for(i = 0; i < n; i++) {
    k = idptr[i] - 1;
    fitptr[i] = fitrptr[k];
    etaptr[i] += fitptr[i];
  }

  /* Evaluate loglik, weights and score vector. */
  peta = map2par(getListElement(family, "map2par"), eta2, rho);
  double pibetaprop = REAL(iwls_eval(VECTOR_ELT(family, ll_ind), response, peta, rho))[0];
  weights = iwls_eval(getListElement(getListElement(family, "weights"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho);
  weightsptr = REAL(weights);
  score = iwls_eval(getListElement(getListElement(family, "score"),
    CHAR(STRING_ELT(id, 0))), response, peta, rho);
  scoreptr = REAL(score);

  xweightsptr[0] = 0;
  xrresptr[0] = 0;

  j = 0;
  for(i = 0; i < n; i++) {
    if(indptr[i] > (j + 1)) {
      for(jj = 0; jj < nc; jj++) {
        XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
      }
      ++j;
      xweightsptr[j] = 0;
      xrresptr[j] = 0;
    }
    k = orderptr[i] - 1;

    zptr[k] = etaptr[k] + scoreptr[k] / weightsptr[k];
    etaptr[k] -= fitptr[k];
    eptr[k] = zptr[k] - etaptr[k];
    xweightsptr[j] += weightsptr[k];
    xrresptr[j] += weightsptr[k] * eptr[k];
  }

  for(jj = 0; jj < nc; jj++) {
    XWptr[jj + nc * j] = Xptr[j + nr * jj] * xweightsptr[j];
  }

  /* Compute X'WX. */
  F77_CALL(dgemm)(transa, transb, &nc, &nc, &nr, &one,
    XWptr, &nc, Xptr, &nr, &zero, XWXptr, &nc);

  /* Add penalty matrix and variance parameter. */
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      for(i = 0; i < nc; i++) {
        for(j = 0; j < nc; j++) {
          XWXptr[i + nc * j] += 1 / tau2ptr[jj] * Sptr[i + nc * j];
        }
      }
    }
  }

  /* Cholesky decompostion of XWX. */
  SEXP L2;
  PROTECT(L2 = duplicate(getListElement(x, "XWX")));
  ++nProtected;
  Lptr = REAL(L2);
  for(j = 0; j < nc; j++) { 	/* Zero the lower triangle. */
	  for(i = j + 1; i < nc; i++) {
      Lptr[i + nc * j] = 0.0;
    }
  }

	F77_CALL(dpotrf)("Upper", &nc, Lptr, &nc, &info);

  /* Compute the inverse precision matrix. */
  SEXP PINV2;
  PROTECT(PINV2 = duplicate(L2));
  ++nProtected;
  PINVptr = REAL(PINV2);

  F77_CALL(dpotri)("Upper", &nc, PINVptr, &nc, &info);
	F77_CALL(dpotrf)("Upper", &nc, PINVLptr, &nc, &info);

  sdiag0 = 0.0;
	for(j = 0; j < nc; j++) {
    sdiag0 += log(pow(Lptr[j + nc * j], 2.0));
	  for(i = j + 1; i < nc; i++) {
		  PINVptr[i + j * nc] = PINVptr[j + i * nc];
    }
  }

  /* Compute mu. */
  F77_CALL(dgemm)(transa2, transb, &nc, &k1, &nr, &one,
    Xptr, &nr, xrresptr, &nr, &zero, mu0ptr, &nc);

  F77_CALL(dgemm)(transa, transb, &nc, &k1, &nc, &one,
    PINVptr, &nc, mu0ptr, &nc, &zero, mu1ptr, &nc);

  /* Log priors. */
  double qbeta = 0.0;
  double p1 = 0.0;
  double a2 = 0.0;
  double b2 = 0.0;
  if(fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      Sptr = REAL(VECTOR_ELT(VECTOR_ELT(x, S_ind), jj));
      gSg = 0.0;
      for(i = 0; i < nc; i++) {
        tsum1 = 0.0;
        tsum2 = 0.0;
        for(j = 0; j < nc; j++) {
          tsum1 += thetaptr[j] * Sptr[j + i * nc];
          tsum2 += (thetaptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
        }
        gSg += tsum1 * thetaptr[i];
        qbeta += tsum2 * (thetaptr[i] - mu1ptr[i]);
      }
      p1 += -log(tau2ptr[jj]) * rankptr[jj] / 2 + -0.5 / tau2ptr[jj] * gSg +
        log(pow(b, a)) - exp(lgamma(a)) + (-a - 1) * log(tau2ptr[jj]) - b / tau2ptr[jj];
      if(fxsp < 1) {
        a2 = rankptr[jj] / 2 + a;
        b2 = 0.5 * gSg + b;
        GetRNGstate();
        tau2ptr[jj] = 1 / rgamma(a2, 1 / b2);
        PutRNGstate();
      }
    }
  } else {
    for(i = 0; i < nc; i++) {
      tsum1 = 0.0;
      p1 += dnorm(thetaptr[i], 0, 1000, 1);
      for(j = 0; j < nc; j++) {
        tsum1 += (thetaptr[j] - mu1ptr[j]) * XWXptr[j + i * nc];
      }
      qbeta += tsum1 * (thetaptr[i] - mu1ptr[i]);
    }
  }
  qbeta = 0.5 * sdiag0 - 0.5 * qbeta;

  SEXP alpha;
  PROTECT(alpha = allocVector(REALSXP, 1));
  ++nProtected;
  REAL(alpha)[0] = (pibetaprop + qbeta + p2) - (pibeta + qbetaprop + p1);

/*Rprintf("pibetaprop %g\n", pibetaprop);*/
/*Rprintf("qbeta %g\n", qbeta);*/
/*Rprintf("p2 %g\n", p2);*/
/*Rprintf("pibeta %g\n", pibeta);*/
/*Rprintf("qbetaprop %g\n", qbetaprop);*/
/*Rprintf("p1 %g\n", p1);*/
/*Rprintf("alpha %g\n", exp(REAL(alpha)[0]));*/

  /* Stuff everything together. */
  SEXP rval;
  PROTECT(rval = allocVector(VECSXP, 2));
  ++nProtected;

  for(j = 0; j < nc; j++) {
    thetaptr[j] = gamma1ptr[j];
  }
  if(fxsp < 1 && fixed < 1) {
    for(jj = 0; jj < ntau2; jj++) {
      thetaptr[nc + jj] = tau2ptr[jj];
    }
  }

  SET_VECTOR_ELT(rval, 0, theta2);
  SET_VECTOR_ELT(rval, 1, alpha);

  SEXP nrval;
  PROTECT(nrval = allocVector(STRSXP, 2));
  ++nProtected;

  SET_STRING_ELT(nrval, 0, mkChar("parameters"));
  SET_STRING_ELT(nrval, 1, mkChar("alpha"));
        
  setAttrib(rval, R_NamesSymbol, nrval);

  UNPROTECT(nProtected);
  return rval;
}


/* Censored normal: left = 0, right = Inf. */
SEXP cnorm_loglik(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, 1));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  int *checkptr = INTEGER(check);

  double ll = 0.0;
  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ll += pnorm5((-1.0 * muptr[i]) / sigmaptr[i], 0.0, 1.0, 1, 1);
    } else {
      ll += dnorm((yptr[i] - muptr[i]) / sigmaptr[i], 0.0, 1.0, 1) - log(sigmaptr[i]);
    }
  }

  REAL(rval)[0] = ll;
  UNPROTECT(1);
  return rval;
}


SEXP cnorm_score_mu(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      rvalptr[i] = -1 * mills / sigmaptr[i];
    } else {
      rvalptr[i] = (yptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0);
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cnorm_score_sigma(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      rvalptr[i] = mills * muptr[i] / sigmaptr[i];
    } else {
      rvalptr[i] = (yptr[i] - muptr[i]) / pow(sigmaptr[i], 2.0) * (yptr[i] - muptr[i]) - 1.0;
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cnorm_weights_mu(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills, d1, d2;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      d1 = -muptr[i] / pow(sigmaptr[i], 2.0);
      d2 = d1 * -muptr[i];
      rvalptr[i] = -1 * (-d1 / sigmaptr[i] * mills - pow(mills, 2.0) / pow(sigmaptr[i], 2.0));
    } else {
      rvalptr[i] = 1 / pow(sigmaptr[i], 2.0);
    }
  }

  UNPROTECT(1);
  return rval;
}


SEXP cnorm_weights_sigma(SEXP y, SEXP mu, SEXP sigma, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  int *checkptr = INTEGER(check);

  double ddist, pdist, mills, d1, d2;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ddist = dnorm(-muptr[i] / sigmaptr[i], 0.0, 1.0, 0) / sigmaptr[i];
      pdist = pnorm5(-muptr[i] / sigmaptr[i], 0.0, 1.0, 1, 0);
      mills = sigmaptr[i] * ddist / pdist;
      d1 = -muptr[i] / pow(sigmaptr[i], 2.0);
      d2 = d1 * -muptr[i];
      rvalptr[i] = -1 * ((-muptr[i] / sigmaptr[i] - (0.0 - muptr[i]) * d2) * mills - pow(muptr[i], 2.0) / pow(sigmaptr[i], 2.0) * pow(mills, 2.0));
    } else {
      rvalptr[i] = 2 / pow(sigmaptr[i], 2.0) * pow(yptr[i] - muptr[i], 2.0);
    }
  }

  UNPROTECT(1);
  return rval;
}


/* Censored normal, left = 0, right - Inf, with power parameter. */
SEXP cnorm_power_loglik(SEXP y, SEXP mu, SEXP sigma, SEXP alpha, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, 1));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double *sigmaptr = REAL(sigma);
  double *alphaptr = REAL(alpha);
  int *checkptr = INTEGER(check);

  double ll = 0.0;
  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      ll += pnorm5(0.0 , muptr[i], sigmaptr[i], 1, 1);
    } else {
      ll += dnorm(pow(yptr[i], 1.0 / alphaptr[i]), muptr[i], sigmaptr[i], 1) -
        log(alphaptr[i]) + (1.0 / alphaptr[i] - 1.0) * log(yptr[i]);
    }
  }

  REAL(rval)[0] = ll;
  UNPROTECT(1);
  return rval;
}


SEXP cnorm_power_score_alpha(SEXP y, SEXP mu, SEXP sigma, SEXP alpha, SEXP check)
{
  SEXP rval;
  PROTECT(rval = allocVector(REALSXP, length(y)));
  int i;
  int n = length(y);
  double *yptr = REAL(y);
  double *muptr = REAL(mu);
  double * sigmaptr = REAL(sigma);
  double *rvalptr = REAL(rval);
  double *alphaptr = REAL(alpha);
  int *checkptr = INTEGER(check);

  double tmp, tmp2, tmp3;

  for(i = 0; i < n; i++) {
    if(checkptr[i]) {
      rvalptr[i] = 0.0;
    } else {
      tmp = exp(-log(alphaptr[i]));
      tmp2 = log(yptr[i]);
      tmp3 = pow(yptr[i], tmp);
      rvalptr[i] = (tmp * tmp3 * tmp2 * (tmp3 - muptr[i])) / pow(sigmaptr[i], 2.0) - tmp * tmp2 - 1.0;
    }
  }

  UNPROTECT(1);
  return rval;
}


/* Fast quantile computation */
void swapd(double *a, double *b)      
{ 
  double temp;
  temp = *a;
  *a = *b;
  *b = temp;
} 
    
void quicksort_body(double *x, int up, int down)
{ 
  int start, end;
  start = up;              
  end = down;                
  while(up < down) {            
    while(x[down] >= x[up] && up < down)      
      down--;              
    if(up != down) {                    
      swapd(&x[up], &x[down]);   
      up++;    
    }
    while(x[up] <= x[down] && up < down)        
      up++;                 
    if(up != down) {                    
      swapd(&x[up], &x[down]);   
      down--;   
    } 
  }       
  if(start < up)   
    quicksort_body(x, start, up - 1); 
  if(end > down)  
    quicksort_body(x, down + 1, end);  
}

void quicksort(int n, double *x)
{ 
  quicksort_body(x, 0, n - 1);    
}


SEXP quick_quantiles(SEXP X, SEXP samples)
{
  int i, j, ii;
  int iter, nr, nc, nProtected = 0;
  SEXP out, TMP, q1, q2, q3, names;
    
  nr = nrows(X);
  nc = ncols(X);
  iter = ncols(samples);
    
  PROTECT(names = allocVector(STRSXP, 3));
  ++nProtected;
        
  PROTECT(out = allocVector(VECSXP, 3));
  ++nProtected;

  PROTECT(TMP = allocVector(REALSXP, iter));
  ++nProtected;
    
  PROTECT(q1 = allocVector(REALSXP, nr));
  ++nProtected;
    
  PROTECT(q2 = allocVector(REALSXP, nr));
  ++nProtected;
    
  PROTECT(q3 = allocVector(REALSXP, nr));
  ++nProtected;
    
  double np11 = iter * 0.025;
  double np12 = iter * 0.5;
  double np13 = iter * 0.975;
    
  int np1 = iter - np11;
  int np2 = iter - np12;
  int np3 = iter - np13;
    
  double *Xptr, *sptr, *tptr, *q1ptr,*q2ptr, *q3ptr;
  Xptr = REAL(X);
  sptr = REAL(samples);
  tptr = REAL(TMP);
  q1ptr = REAL(q1);
  q2ptr = REAL(q2);
  q3ptr = REAL(q3);

  double tmp = 0.0;
    
  for(i = 0; i < nr; ++i) {
    for(ii = 0; ii < iter; ++ii) {
      tmp = 0.0;
      for(j = 0; j < nc; ++j) {
        tmp += Xptr[i + j * nr] * sptr[j + nc * ii];
      }
      tptr[ii] = tmp;
    }

    quicksort(iter, tptr);
              
    if((np11 - floor(np11)) == 0.0) {
      q1ptr[i] = (tptr[np1 - 1] + tptr[np1]) / 2.0;
    } else {
      q1ptr[i] = tptr[np1 - 1];
    }
    if((np12 - floor(np12)) == 0.0) { 
      q2ptr[i] = (tptr[np2 - 1] + tptr[np2]) / 2.0;
    } else {
      q2ptr[i] = tptr[np2 - 1];
    }
    if((np13 - floor(np13)) == 0.0) { 
      q3ptr[i] = (tptr[np3 - 1] + tptr[np3]) / 2.0;
    } else {
      q3ptr[i] = tptr[np3 - 1];
    }
  }
        
  SET_VECTOR_ELT(out, 0, q1);
  SET_VECTOR_ELT(out, 1, q2);
  SET_VECTOR_ELT(out, 2, q3);
    
  SET_STRING_ELT(names, 0, mkChar("lo"));
  SET_STRING_ELT(names, 1, mkChar("med"));
  SET_STRING_ELT(names, 2, mkChar("up"));
    
  setAttrib(out, R_NamesSymbol, names);
    
  UNPROTECT(nProtected);
  return out;
}

