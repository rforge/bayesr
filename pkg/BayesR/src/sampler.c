#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Applic.h> /* for dgemm */

#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>
#include <R_ext/RS.h> 
#include <R_ext/BLAS.h> 
#include <R_ext/Lapack.h> 
#include <R_ext/Linpack.h> 

#define getDims(A) INTEGER(coerceVector(getAttrib(A, R_DimSymbol), INTSXP))
    
     
SEXP pvec(SEXP vec)
    {
    int i,n = length(vec);
    double *vecptr = REAL(vec);
    for(i=0;i<n;i++,vecptr++)
        Rprintf(" %g",*vecptr);
    Rprintf("\n");
    return R_NilValue;
    }
    
SEXP pmat(SEXP mat)
    {
    int i,j,n = nrows(mat),k=ncols(mat);
    Rprintf("   ");
    for(j=0;j<k;++j)
        Rprintf("[%d] ",j);
    Rprintf("\n");
    for(i=0;i<n;++i)
        {
        Rprintf("[%d]",i);
        for(j=0;j<k;++j)
            Rprintf(" %g",REAL(mat)[i + j*n]);
        Rprintf("\n");
        }
    Rprintf("\n");
    return R_NilValue;
    }
    
/* oldstuff...       
SEXP vadd(SEXP A, SEXP B) // elementwise vector addition
    {
    int i,n;
    
    n = length(A);
    
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    
    double *ap,*bp,*op;
    ap = REAL(A);
    bp = REAL(B);
    op = REAL(out);
    
    for(i=0;i<n;++i)
        op[i] = ap[i] + bp[i];
        
    UNPROTECT(1);
    return out;
    }
    
    
void vaddc(double *a, double *b, int n) // elementwise vector addition
    {
    int i;
    
    for(i=0;i<n;++i)
        a[i] += b[i];
    }


SEXP vmin(SEXP A, SEXP B) // elementwise vector substraction
    {
    int i,n;
    
    n = length(A);
    
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    
    double *ap,*bp,*op;
    ap = REAL(A);
    bp = REAL(B);
    op = REAL(out);
    
    for(i=0;i<n;++i)
        op[i] = ap[i] - bp[i];
        
    UNPROTECT(1);
    return out;
    }
    
    
void vminc(double *a, double *b, int n) // elementwise vector substraction
    {
    int i;
    
    for(i=0;i<n;++i)
        a[i] -= b[i];
    }

    
    
SEXP vmul(SEXP A, SEXP B) // elementwise vector multiplication
    {
    int i,n;
    
    n = length(A);
    
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    
    double *ap,*bp,*op;
    ap = REAL(A);
    bp = REAL(B);
    op = REAL(out);
    
    for(i=0;i<n;++i)
        op[i] = ap[i] * bp[i];
        
    UNPROTECT(1);
    return out;
    }
    
    
void vmulc(double *a, double *b, int n) // elementwise vector multiplication
    {
    int i; 
    
    for(i=0;i<n;++i)
        a[i] *= b[i];
    }

    
    
SEXP vmuls(SEXP A, double b) // elementwise vector multiplication with scalar A*b
    {
    int i,n;
    
    n = length(A);
    
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    
    double *ap,*op;
    ap = REAL(A);
    op = REAL(out);
    
    for(i=0;i<n;++i)
        op[i] = ap[i] * b;
        
    UNPROTECT(1);
    return out;
    }
    
    
void vmulsc(double *a, double b, int n) // elementwise vector multiplication with scalar
    {
    int i;
    
    for(i=0;i<n;++i)
        a[i] *= b;
    }
    
    
SEXP vmulh(double a, SEXP B) // special help multiplication
    {
    int i,n;   
    
    n = length(B);
    
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    
    double *bp,*op;
    bp = REAL(B);
    op = REAL(out);
    
    for(i=0;i<n;++i)
        op[i] = 1/(1 + a * bp[i]);
        
    UNPROTECT(1);
    return out;
    }
    
    
void vmulhc(double a, double *b, int n) // special help multiplication
    {
    int i;
    
    for(i=0;i<n;++i)
        b[i] = 1/(1 + a * b[i]);
    }
        

SEXP mac(SEXP mat, int k) // get column k of matrix mat
    {
    int i,n;
    n = nrows(mat);
      
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    
    double *outp,*mp;
    outp = REAL(out);
    mp = REAL(mat);
    
    for(i=0;i<n;++i)
        {
        outp[i] = mp[k*n + i];
        }
        
    UNPROTECT(1);
    return out;
    }
    
    
void macc(double *a, double *b, int k, int n)
    {
    int i;
    
    for(i=0;i<n;++i)
        b[i] = a[k*n + i];
    }

    
void smac(SEXP mat, SEXP what, int k) // set column k of matrix mat
    {
    int i,n;    
    n = nrows(mat);
    
    double *wp,*mp;
    wp = REAL(what);
    mp = REAL(mat);
    
    for(i=0;i<n;++i)
        {
        mp[k*n + i] = wp[i];
        }     
    }
    
    
void smacc(double *mat, double *what, int k, int n)
    {
    int i;
    for(i=0;i<n;++i)
        mat[k*n + i] = what[i];
    }

    
    
SEXP gs(SEXP list, int what) // get the level 2 list element
    {
    SEXP elmt = R_NilValue;
  
    elmt = VECTOR_ELT(list, what);
    elmt = VECTOR_ELT(elmt, 1);

    return elmt;
    }
    
    
SEXP mprod(SEXP x, SEXP y) // simple matrix %*% matrix product
    {
    int nrx, ncx, nry, ncy;
    int *xdims, *ydims;
    double *ansptr, *xptr, *yptr;

    xdims = getDims(x);  
    ydims = getDims(y);
    nrx = xdims[0];      
    ncx = xdims[1];
    nry = ydims[0];      
    ncy = ydims[1];

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
    

SEXP cvp(SEXP x, SEXP y) //crossproduct matrix with vector
    {
    int nrx, ncx, nry, ncy;
    int *xdims;
    double *ansptr, *xptr, *yptr;

    xdims = getDims(x);  
    nrx = xdims[0];      
    ncx = xdims[1];
    nry = nrx;      
    ncy = 1;

    //PROTECT(x = coerceVector(x, REALSXP));
    //PROTECT(y = coerceVector(y, REALSXP));
    xptr = REAL(x);  
    yptr = REAL(y);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, ncx));
    ansptr = REAL(ans);
  
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    
    if (nrx>0 && ncx>0 && nry>0 && ncy>0) 
        {
        F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
                        xptr, &nrx, yptr, &nry, &zero, ansptr, &ncx);
        }
    else
        {
        int i;
        for(i=0;i<ncx*ncy;i++) 
            ansptr[i] = 0;
        }
    UNPROTECT(1);
    return(ans);
    }   
    
    
SEXP vrnorm(SEXP mu, SEXP s2) // draw random numbers for a vector from the normal distribution
    {
    //int nProtected = 0;
    //PROTECT(mu = coerceVector(mu, REALSXP));
    //++nProtected;
    //PROTECT(s2 = coerceVector(s2, REALSXP));
    //++nProtected;
    
    int i,n;
    n = length(mu);
      
    SEXP out;
    PROTECT(out = allocVector(REALSXP, n));
    //++nProtected;
    
    double *muptr,*s2ptr,*optr,tmp;
    muptr = REAL(mu);
    s2ptr = REAL(s2);
    optr = REAL(out);
    
    GetRNGstate();
    for(i=0;i<n;++i)
        {
        tmp = sqrt(s2ptr[i]);
        optr[i] = rnorm(muptr[i], tmp);
        }
    PutRNGstate();
     
    UNPROTECT(1);
    return out;   
    } 
    
    
double sumit(SEXP x) // conventional sum
    {
    int i,n;
    n = length(x);
    
    double sum = 0.0,*xptr;
    xptr = REAL(x);
    
    for(i=0;i<n;++i)
        sum += xptr[i];
        
    return sum;
    }
    

void setmatelmt(SEXP mat, double what, int row, int col) // set element i,j of matrix mat
    {
    int n = nrows(mat);
    
    REAL(mat)[col*n + row] = what;
    }
    
    
double ssum(SEXP x) // calculate sum of squared elements of vector
    {
    int i,n = length(x);
    double sum = 0,*xptr;
    xptr = REAL(x);
    
    for(i=0;i<n;++i)
        {
        sum += pow(xptr[i],2);
        }
        
    return sum;
    }
    
    
SEXP ngamma(SEXP a, SEXP b)
    {
    SEXP out;
    PROTECT(out = allocVector(REALSXP,1));
    
    REAL(out)[0] = rgamma(REAL(a)[0],REAL(b)[0]);
    
    UNPROTECT(1);
    return out;
    }
    
  
double ngma(double a, double b)
    {
    double out;    
    out = rgamma(a,b);
    return out;
    }
*/  
    
void cvp2(SEXP x, SEXP y, SEXP ans) //crossproduct matrix with vector
    {
    int nrx, ncx, nry, ncy;
    int *xdims;
    double *ansptr, *xptr, *yptr;

    xdims = getDims(x);  
    nrx = xdims[0];      
    ncx = xdims[1];
    nry = nrx;      
    ncy = 1;

    xptr = REAL(x);  
    yptr = REAL(y);

    ansptr = REAL(ans);
  
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    
    if (nrx>0 && ncx>0 && nry>0 && ncy>0) 
        {
        F77_CALL(dgemm)(transa, transb, &ncx, &ncy, &nrx, &one,
                        xptr, &nrx, yptr, &nry, &zero, ansptr, &ncx);
        }
    else
        {
        int i;
        for(i=0;i<ncx*ncy;i++) 
            ansptr[i] = 0;
        }
    } 
    
    
void cvp1(SEXP x, SEXP y, SEXP ans) //product matrix with vector
    {
    int nrx, ncx, nry, ncy;
    int *xdims;
    double *ansptr, *xptr, *yptr;

    xdims = getDims(x);  
    nrx = xdims[0];      
    ncx = xdims[1];
    nry = ncx;      
    ncy = 1;

    xptr = REAL(x);  
    yptr = REAL(y);

    ansptr = REAL(ans);
  
    char *transa = "N", *transb = "N";
    double one = 1.0, zero = 0.0;
    
    if (nrx>0 && ncx>0 && nry>0 && ncy>0) 
        {                       
        F77_CALL(dgemm)(transa, transb, &nrx, &ncy, &ncx, &one,
                        xptr, &nrx, yptr, &nry, &zero, ansptr, &nrx);
        }
    else
        {
        int i;
        for(i=0;i<ncx*ncy;i++) 
            ansptr[i] = 0;
        }
    } 


double sums(SEXP beta, SEXP K)
    {
    int nrb, ncb, nrk, nck, i;
    int *kdims;
    double *ansptr, *bptr, *kptr, out;

    kdims = getDims(K);  
    nrk = kdims[0];      
    nck = kdims[1];
    nrb = nrk;      
    ncb = 1;

    bptr = REAL(beta);  
    kptr = REAL(K);

    SEXP ans;
    PROTECT(ans = allocVector(REALSXP, nrb));
    ansptr = REAL(ans);
  
    char *transa = "T", *transb = "N";
    double one = 1.0, zero = 0.0;
    
    if (nrb>0 && ncb>0 && nrk>0 && nck>0) 
        {
        F77_CALL(dgemm)(transa, transb, &ncb, &nck, &nrb, &one,
                        bptr, &nrb, kptr, &nrk, &zero, ansptr, &ncb);
        }
    else
        {
        for(i=0;i<ncb*nck;i++) 
            ansptr[i] = 0;
        }
    for(i=0;i<ncb*nck;i++)
	out += ansptr[i]*bptr[i];

    UNPROTECT(1);
    return out;
    }
    
    
void gcenter(SEXP x) // vector centering
    {
    int ii, N;
    double mymean, mysum = 0, *xptr;
    
    //PROTECT(x = coerceVector(x, REALSXP));
	
    xptr = REAL(x);
    N = length(x);

    for(ii = 0; ii < N; ++ii)
        {
        mysum += xptr[ii];
        }

    mymean = mysum/N;
	
    for(ii = 0; ii < N; ++ii)
        {
        xptr[ii] = xptr[ii] - mymean;
        }
    }
    
   
int thincheck(SEXP iterthin, int which) // check for storing samples
    {
    int n,j,check = 0;
    
    n = length(iterthin);
    int *tmp = INTEGER(iterthin); 
    
    for(j=0;j<n;++j)
        {       
        if(which == (tmp[j]))
            {
            check = 1;
            break;
            }
        }

    return check;
    }
   
    
double dev(SEXP response, SEXP eta, double sigma2, int n)
    {
    int i;
    double tmp1,tmp2,*rptr,*eptr;
    rptr = REAL(response);
    eptr = REAL(eta);
    
    tmp2 = 0;
    tmp1 = 0;
        
    for(i=0;i<n;++i)
        {
        tmp1 = (rptr[i] - eptr[i]);
        tmp1 = tmp1*tmp1;
        tmp2 += log(2*M_PI*sigma2) + 1/sigma2 * tmp1;
        }
        
    return tmp2;
    }
    
    
double tnorm(double low, double high, double mean, double sigma)
    {
    GetRNGstate();
    double p = pnorm(low,mean,sigma,1,0) + unif_rand()*(pnorm(high,mean,sigma,1,0) - pnorm(low,mean,sigma,1,0));
    PutRNGstate();
    double out = qnorm(p,mean,sigma,1,0);
    return out;
    }
    
    
void vtnorm(SEXP low, SEXP high, SEXP mean, double sigma, SEXP out, int n)
    {
    int i;
    double *lowptr,*highptr,*meanptr,*outptr;
    lowptr = REAL(low);
    highptr = REAL(high);
    meanptr = REAL(mean);
    outptr = REAL(out);
    
    for(i=0;i<n;++i)
        outptr[i] = tnorm(lowptr[i],highptr[i],meanptr[i],sigma);
    }


double Phi2(double x)
  	{
  	double xt;
  	if (x < 0)
   		xt = -x;
  	else
    		xt = x;
  	double t=xt*xt;
  	double pol = 1+0.196854*xt+0.115194*t;
  	t*=xt;
  	pol+=0.000344*t;
  	t*=xt;
  	pol+=0.019527*t;

  	pol = 1.0/pow(pol,4);

  	if (x < 0)
    		return 0.5*pol;
  	else
    		return 1.0-0.5*pol;
  	}


double invPhi2(double p)
	{
    	double pt;
    	if(p < 0.5)
      		pt = p;
    	else
      		pt = 1-p;

    	if(pt < 1.0e-100)
      		pt = 1.0e-100;

    	double t = sqrt(-2.0*log(pt));
    	double t2 = t*t;
    	double t3 = t2*t;
    	double r1 = 2.515517+0.802853*t+0.010328*t2;
    	double r2 = 1+1.432788*t+0.189269*t2+0.001308*t3;

    	if(p<0.5)
      		return -t+r1/r2;
    	else
      		return t-r1/r2;
  	}


double trunc_normal2(double a, double b, double mu, double s)
	{
  	double at = Phi2((a-mu)/s);
 	double bt = Phi2((b-mu)/s);
	GetRNGstate();
  	double u = at+(bt-at)*unif_rand();
	PutRNGstate();
  	double r = mu+s*invPhi2(u);
  	if(r < a)
    		r = a+0.00000001;
  	if(r > b)
    		r = b-0.00000001;
  	return r;
  	}


SEXP rtvnorm(SEXP n, SEXP mean, SEXP sd, SEXP low, SEXP high)
	{
	int i;
    	double *lowptr,*highptr,*meanptr,*outptr,*sdptr;
    	lowptr = REAL(low);
    	highptr = REAL(high);
    	meanptr = REAL(mean);
	sdptr = REAL(sd);

	SEXP out;
    	PROTECT(out = allocVector(REALSXP, INTEGER(n)[0]));
    	outptr = REAL(out);

    	for(i=0;i<INTEGER(n)[0];++i)
        	outptr[i] = tnorm(lowptr[i],highptr[i],meanptr[i],sdptr[i]);

    	UNPROTECT(1);
    	return out;
	}


SEXP rtvnorm2(SEXP n, SEXP mean, SEXP sd, SEXP low, SEXP high)
	{
	int i;
    	double *lowptr,*highptr,*meanptr,*outptr,*sdptr;
    	lowptr = REAL(low);
    	highptr = REAL(high);
    	meanptr = REAL(mean);
	sdptr = REAL(sd);

	SEXP out;
    	PROTECT(out = allocVector(REALSXP, INTEGER(n)[0]));
    	outptr = REAL(out);

    	for(i=0;i<INTEGER(n)[0];++i)
        	outptr[i] = trunc_normal2(lowptr[i],highptr[i],meanptr[i],sdptr[i]);

    	UNPROTECT(1);
    	return out;
	}


double rand_expo(double lambda)
  	{
	GetRNGstate();
  	return (- 1/lambda)*log(unif_rand());
	PutRNGstate();
  	}


double rand_gamma(double a,double b)
	{
	GetRNGstate();
  	if (a > 1)
	 	{
	 	double h1 = a-1;         // h1 entspricht b in Devroye (1986)
	 	double h2 = 3*a-0.75;    // h2 entspricht c in Devroye (1986)
	 	double U,V,W,Y,X,Z;
	 	int accept = 0;
	 	do
			{
			U = unif_rand();
			V = unif_rand();
			W = U*(1-U);
			Y = sqrt(h2/W)*(U-0.5);
			X = h1 + Y;
			if (X > 0)
				{
				Z = 64*W*W*W*V*V;
				if ( Z <= (1 - (2*Y*Y)/X) )
			  		accept = 1;
				else
			  		{
			  		if ( ((X/h1) > 0) &&  ( log(Z) <= ( 2*(h1*log(X/h1) - Y) ) ) )
				 		accept = 1;
			  		}
				}
			}
	 	while (accept == 0);
	 		return X/b;
	 	}
  	else
	 	{
	 	if (a == 1)
			return rand_expo(b);
	 	else
			{
			double X = rand_gamma(a+1,1)*pow(unif_rand(),1/a);
			return X/b;
			}
	 	}
	PutRNGstate();
  	}


SEXP rgamma2(SEXP n, SEXP a, SEXP b)
	{
	int i;
	double *outptr;

	SEXP out;
    	PROTECT(out = allocVector(REALSXP, INTEGER(n)[0]));
    	outptr = REAL(out);

	for(i=0;i<INTEGER(n)[0];++i)
		{
		GetRNGstate();
		outptr[i] = rand_gamma(REAL(a)[0],REAL(b)[0]);
		PutRNGstate();
		}

	UNPROTECT(1);
	return out;
	}


void sigmat(SEXP ZZ, SEXP K, double lambda, SEXP out, int dim)
	{
	int i,j;
	double *ZZptr=REAL(ZZ),*Kptr=REAL(K),*outptr=REAL(out);
	for(i=0;i<dim;++i)
		for(j=0;j<dim;++j)
			outptr[i + dim*j] = ZZptr[i + dim*j] + lambda*Kptr[i + dim*j];
	}


void monysamp(int maxit, 
              SEXP ZZ, 
              SEXP K, 
              double lambda, 
              int dP,
              SEXP basis, 
              SEXP e, 
              SEXP mus, 
              SEXP draws, 
              double sumtmp, 
              double sigma2, 
              int check)
	{
	int mcount = 0,ii,dd;
	double mcsmu;

	SEXP sigmats;
	PROTECT(sigmats = allocVector(REALSXP,dP*dP));

	sigmat(ZZ,K,lambda,sigmats,dP);
	cvp2(basis,e,mus);	

	double *sigmatsptr = REAL(sigmats);
	double *musptr = REAL(mus);
	double *drawssptr = REAL(draws);

	while(mcount < maxit)
		{
		for(ii=0;ii<dP;++ii)
			{
			mcsmu = 0;

			if(ii>0)
				for(dd=0;dd<ii;++dd)
					mcsmu += drawssptr[dd]*sigmatsptr[ii + dP*dd];
			if(ii<(dP - 1))
				for(dd=(ii+1);dd<dP;++dd)
					mcsmu += drawssptr[dd]*sigmatsptr[ii + dP*dd];
				
			mcsmu = (musptr[ii] - mcsmu)/sigmatsptr[ii + dP*ii];

			if(check == 1)
				{
				if(ii==0)
					drawssptr[ii] = trunc_normal2(-1000,drawssptr[1],mcsmu,sqrt(sigma2/sigmatsptr[ii + dP*ii]));
				if(ii==(dP-1))
					drawssptr[ii] = trunc_normal2(drawssptr[ii-1],1000,mcsmu,sqrt(sigma2/sigmatsptr[ii + dP*ii]));
				if(ii>0 && ii < (dP-1)) 
					drawssptr[ii] = trunc_normal2(drawssptr[ii-1],drawssptr[ii+1],mcsmu,sqrt(sigma2/sigmatsptr[ii + dP*ii]));
				}
			else
				{
				if(ii==0)
					drawssptr[ii] = trunc_normal2(drawssptr[1],1000,mcsmu,sqrt(sigma2/sigmatsptr[ii + dP*ii]));
				if(ii==(dP-1))
					drawssptr[ii] = trunc_normal2(-1000,drawssptr[ii-1],mcsmu,sqrt(sigma2/sigmatsptr[ii + dP*ii]));
				if(ii>0 && ii < (dP-1)) 
					drawssptr[ii] = trunc_normal2(drawssptr[ii+1],drawssptr[ii-1],mcsmu,sqrt(sigma2/sigmatsptr[ii + dP*ii]));
				}
			}
		mcount++;
		}

	sumtmp = sums(draws, K);
	UNPROTECT(1);
	}


SEXP rtr2(SEXP help, SEXP coefs)
	{
	int i,np = length(coefs);
	double u = 0.0,sc = 0.0;

	SEXP out;
	PROTECT(out = allocVector(REALSXP,np));

	double *h = REAL(help);
	double *c = REAL(coefs);
	double *o = REAL(out);
	
	for(i=0;i<np;++i)
		{
		u += h[i];
		sc += c[i];
		}
	for(i=0;i<np;++i)
		o[i] = c[i] - (1/u * h[i])*sc;

	UNPROTECT(1);
	return out;
	}


void rtr(double *s, double *c, double la, int np, double sum)
	{
	int i;
	double u = 0.0,sc = 0.0;
	sum = 0.0;

	for(i=0;i<np;++i)
		{
		u += 1/(1+la*s[i]);
		sc += c[i];
		}
	for(i=0;i<np;++i)
		{
		c[i] = c[i] - (1/u * (1/(1+la*s[i])))*sc;
		sum += (c[i]*c[i])*s[i];
		}
	}
    
    
// random sampling 3rd stage and more
void rsampling(SEXP ins, 
                SEXP response, 
                SEXP eta, 
                SEXP e,
                int n, 
                double sig2, 
                double hypbt2, 
                double hypb, 
                double ce,
                int thinit,
                int itrthin,
                int maxit,
		int all)
    {
    int j,jj,ii,r;
    
    r = INTEGER(VECTOR_ELT(ins, 44))[0];
    
    double lambd=0,helpme=0,mur=0,sigmaR=0,drawtmp=0,tmp=0,sumtmp=0,mymean=0,s2tmp=0;

    double *etaRptr,*responseptr,*fRptr,*eptr,*sigma2Rptr,*cRLptr;
    double *S2Rptr,*cRMSptr,*etaptr,*fptr,*stetaRptr;
    etaRptr = 0;
    stetaRptr = 0;
    responseptr = REAL(response); 
    fRptr = REAL(VECTOR_ELT(ins, 0));
    eptr = REAL(e);
    sigma2Rptr = REAL(VECTOR_ELT(ins, 2));
    cRLptr = REAL(VECTOR_ELT(ins, 32));
    S2Rptr = REAL(VECTOR_ELT(ins, 29));
    cRMSptr = REAL(VECTOR_ELT(ins, 17));
    cRLptr = REAL(VECTOR_ELT(ins, 32));
    etaptr = REAL(eta);
    fptr = REAL(VECTOR_ELT(ins, 68));
    
    int *BRptr,*RRiptr,*RRsptr,*RRlptr,*dPRlptr,*dPRptr,*hspecptr,*centerptr;
    int *newRptr,*rindptr,*bcheckptr,*rccptr;
    BRptr = INTEGER(VECTOR_ELT(ins, 9));
    RRiptr = INTEGER(VECTOR_ELT(ins, 48));
    RRsptr = INTEGER(VECTOR_ELT(ins, 46));
    RRlptr = INTEGER(VECTOR_ELT(ins, 47));
    dPRlptr = INTEGER(VECTOR_ELT(ins, 50));
    hspecptr = INTEGER(VECTOR_ELT(ins, 42)); 
    bcheckptr = INTEGER(VECTOR_ELT(ins, 72));
    rccptr = INTEGER(VECTOR_ELT(ins, 74));
    
    double *metaRptr,*sptr,*cRptr,*cRMptr,*tau2Rptr;
    double *hyperaTau2newRptr,*s1ptr,*frptr;
    double *DRAWERSptr,*DRAWERSMptr,*nilorptr;
    double *DRAWERS1ptr,*DRAWERSM1ptr,*nilor1ptr;
    double *DRAWERS1Lptr,*DRAWERSM1Lptr,*nilor1Lptr;
    double *drawTildeGammaRptr,*drawTildeGammaRMptr;
    double *frlptr,*cRMLptr;
    double *drawTZRlptr,*drawTZRlMptr,*cRLLptr,*drawBptr,*drawBMptr;
    double *drawBetaRptr,*haRnewptr,*cRSptr,*cRMSSptr;
    double *fRtmpptr,*ZZptr;
    
    haRnewptr = REAL(VECTOR_ELT(ins, 25));
    cRptr = REAL(VECTOR_ELT(ins, 31));
    
    for(j=0;j<r;++j)
        {  
        sptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 7), j), 1));
        DRAWERSptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 59), j));
        DRAWERSMptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 60), j));
        drawBetaRptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 10), j));
        nilorptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 61), j));
        cRMptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 6), j));
        rindptr = INTEGER(VECTOR_ELT(VECTOR_ELT(ins, 40), j));
        ZZptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 3), j));
  
        for(ii=0;ii<n;++ii)
            {
            etaptr[ii] = etaptr[ii] - fRptr[n*j + ii];
            eptr[ii] = responseptr[ii] - etaptr[ii];
            }

        cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 7), j), 0), e, VECTOR_ELT(VECTOR_ELT(ins, 61), j));
    
        if(hspecptr[j] > 0)
            {
            etaRptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 37), j));
            metaRptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 4), j));

            s2tmp = 0.0;
                    
            if(LOGICAL(VECTOR_ELT(ins, 41))[j])
                {
                lambd = sig2/sigma2Rptr[j];
                for(ii=0;ii<BRptr[j];++ii)
                    {
                    helpme = 1/(ZZptr[ii] + lambd*sptr[ii]);
                    mur = helpme*(nilorptr[ii] + lambd*etaRptr[ii]);
                    sigmaR = helpme * sig2;
                    sigmaR = sqrt(sigmaR);
                    GetRNGstate();
                    drawtmp = rnorm(mur,sigmaR);
                    PutRNGstate();
	            tmp = drawtmp - etaRptr[ii];
	            s2tmp += tmp*tmp;
                    DRAWERSptr[ii] = drawtmp;
                    DRAWERSMptr[ii] = mur;
                    }
                }
            else
                {
                lambd = sig2/sigma2Rptr[j];
                for(ii=0;ii<BRptr[j];++ii)
                    {
                    tmp = lambd*sptr[ii];
                    helpme = 1/(ZZptr[ii] + tmp);
                    mur = helpme*(nilorptr[ii] + tmp*etaRptr[ii]);
                    sigmaR = helpme * sig2;
                    sigmaR = sqrt(sigmaR);
                    GetRNGstate();
                    drawtmp = rnorm(mur,sigmaR);
                    PutRNGstate();
	            tmp = drawtmp - etaRptr[ii];
	            s2tmp += tmp*tmp;
                    DRAWERSptr[ii] = drawtmp;
                    DRAWERSMptr[ii] = mur;
                    }
                }

/*            s2tmp = 0.0;*/
/*            for(ii=0;ii<BRptr[j];++ii)*/
/*                {*/
/*                tmp = drawBetaRptr[ii] - etaRptr[ii];*/
/*                s2tmp += tmp*tmp;*/
/*                }  */
                        
            // retransform parameters    
            cvp1(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 7), j), 2), VECTOR_ELT(VECTOR_ELT(ins, 59), j), VECTOR_ELT(VECTOR_ELT(ins, 10), j));   

                 

            if(RRiptr[j] > 0)
                {
                rsampling(VECTOR_ELT(VECTOR_ELT(ins, 38), j),VECTOR_ELT(VECTOR_ELT(ins, 10), j),
			  VECTOR_ELT(VECTOR_ELT(ins, 37), j),VECTOR_ELT(VECTOR_ELT(ins, 71), j),
                          BRptr[j],sigma2Rptr[j],hypbt2,hypb,ce,thinit,itrthin,maxit,all);
                }
            if(RRsptr[j] > 0)
                {
                tau2Rptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 13), j));
                dPRptr = INTEGER(VECTOR_ELT(VECTOR_ELT(ins, 15), j));
                hyperaTau2newRptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 22), j));
                centerptr = INTEGER(VECTOR_ELT(VECTOR_ELT(ins, 11), j));
                cRSptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 19), j));
                        
                for(jj=0;jj<RRsptr[j];++jj)
                    {
                    s1ptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 1));
                    frptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 8), j), jj));
                    DRAWERS1ptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 62), j), jj));
                    DRAWERSM1ptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 63), j), jj));
                    nilor1ptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 64), j), jj));  
                    newRptr = INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 70), j), jj));
                    fRtmpptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 69), j), jj));
                            
                    for(ii=0;ii<BRptr[j];++ii)
                        {
                        etaRptr[ii] -= frptr[ii];
                        nilorptr[ii] = drawBetaRptr[ii] - etaRptr[ii];
                        }
                                
                    if(LOGICAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 45), j), jj), 0))[0])
                        {
                        lambd = sigma2Rptr[j]/tau2Rptr[jj];
                        }
                    else
                        {
                        lambd = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 45), j), jj), 1))[0];
                        }

		    if(INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 3))[0] > 0)
	            	{
			monysamp(maxit, 
              			 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 4), 
              			 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 1), 
              			 lambd, 
              			 dPRptr[jj],
              			 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 0), 
              			 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 64), j), jj), 
              			 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 63), j), jj), 
              			 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 62), j), jj), 
              		  	 sumtmp, 
              			 sigma2Rptr[j], 
              			 INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 3))[0]);
			}
		    else
			{                     
                    	cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 14), j), jj), 0), 
                             VECTOR_ELT(VECTOR_ELT(ins, 61), j), VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 64), j), jj));
          
                    	sumtmp = 0.0;                    
                    	for(ii=0;ii<dPRptr[jj];++ii)
                        	{
                        	helpme = 1/(1 + lambd*s1ptr[ii]);
                        	mur = helpme * nilor1ptr[ii];
                        	sigmaR = helpme * sigma2Rptr[j];
                        	sigmaR = sqrt(sigmaR);
                        	GetRNGstate();
                        	drawtmp = rnorm(mur,sigmaR);
                        	PutRNGstate();
                        	DRAWERS1ptr[ii] = drawtmp;
                        	DRAWERSM1ptr[ii] = mur;
                        	sumtmp += (drawtmp*drawtmp)*s1ptr[ii];
                        	}
			}
                               
/*		    if(centerptr[jj] > 0)*/
/*			rtr(s1ptr,DRAWERS1ptr,lambd,dPRptr[jj],sumtmp);*/

                    // sampling smooth variance
                    hypb = hypbt2 + sumtmp/2;
                    tau2Rptr[jj] = 1/rand_gamma(hyperaTau2newRptr[jj],hypb);
                            
                    cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 16), j), jj), VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 62), j), jj), 
                         VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 69), j), jj));
   
                    if(centerptr[jj] > 0)
                        {
                        gcenter(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 69), j), jj));
                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            frptr[ii] = fRtmpptr[newRptr[ii]];
                            etaRptr[ii] += frptr[ii];
                            metaRptr[ii] += etaRptr[ii];
                            }
                        }
                    else
                        {
                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            frptr[ii] = fRtmpptr[newRptr[ii]];
                            etaRptr[ii] += frptr[ii];
                            metaRptr[ii] += etaRptr[ii];
                            }                                
                        }

/*                    for(ii=0;ii<BRptr[j];++ii)*/
/*                    	{*/
/*                        frptr[ii] = fRtmpptr[newRptr[ii]];*/
/*                        etaRptr[ii] += frptr[ii];*/
/*                        metaRptr[ii] += etaRptr[ii];*/
/*                        }    */
                    
                    ce -= cRSptr[jj];
                    cRSptr[jj] = 0;
                    cRMSSptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 18), j), jj), 0));
                    if(LOGICAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 18), j), jj), 1))[0])
                        {
                        cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 18), j), jj), 2),
                             VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 62), j), jj),VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 64), j), jj));
                        for(ii=0;ii<dPRptr[jj];++ii)
                            cRSptr[jj] += cRMSSptr[ii]*nilor1ptr[ii];
                        }
                    else
                        {
                        for(ii=0;ii<dPRptr[jj];++ii)
                            cRSptr[jj] += cRMSSptr[ii]*DRAWERS1ptr[ii];
                        }
                    ce += cRSptr[jj];
                            
                    if(thinit == 1)
                        {
                        REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 26), j), jj))[itrthin] = tau2Rptr[jj];
			if(all > 1)
                        	REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 20), j), jj))[itrthin] = ce - cRSptr[jj];
			else
				REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 20), j), jj))[itrthin] = ce;
                                
                        drawTildeGammaRptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 27), j), jj));
                        drawTildeGammaRMptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 57), j), jj));
                                
                        for(ii=0;ii<dPRptr[jj];++ii)
                            {
                            drawTildeGammaRptr[ii+itrthin*dPRptr[jj]] = DRAWERS1ptr[ii];
                            drawTildeGammaRMptr[ii+itrthin*dPRptr[jj]] = DRAWERSM1ptr[ii];
                            }
                        }
                    }
                }
            if(RRlptr[j] > 0)
                {
                frlptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 55), j));
                DRAWERS1Lptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 65), j));
                DRAWERSM1Lptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 66), j));
                cRMLptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 33), j), 1));
                nilor1Lptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 67), j));
                        
                for(ii=0;ii<BRptr[j];++ii)
                    {
                    etaRptr[ii] -= frlptr[ii];
                    nilorptr[ii] = drawBetaRptr[ii] - etaRptr[ii];
                    }
                cvp2(VECTOR_ELT(VECTOR_ELT(ins, 52), j), VECTOR_ELT(VECTOR_ELT(ins, 61), j), VECTOR_ELT(VECTOR_ELT(ins, 66), j));
                        
                sigmaR = sqrt(sigma2Rptr[j]);
                for(ii=0;ii<dPRlptr[j];++ii)
                    {
	            GetRNGstate();
                    DRAWERS1Lptr[ii] = rnorm(DRAWERSM1Lptr[ii], sigmaR);
                    PutRNGstate();
                    }
                cvp1(VECTOR_ELT(VECTOR_ELT(ins, 52), j),VECTOR_ELT(VECTOR_ELT(ins, 65), j),VECTOR_ELT(VECTOR_ELT(ins, 55), j));

                for(ii=0;ii<BRptr[j];++ii)
                    {
                    etaRptr[ii] += frlptr[ii];
                    metaRptr[ii] += etaRptr[ii];
                    }

                ce -= cRLptr[j];
                cvp1(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 33), j), 0),VECTOR_ELT(VECTOR_ELT(ins, 65), j),VECTOR_ELT(VECTOR_ELT(ins, 67), j));
                tmp = 0;
                for(ii=0;ii<dPRlptr[j];++ii)
                    {
                    nilor1Lptr[ii] *= cRMLptr[ii];
                    tmp += nilor1Lptr[ii];
                    }
                cRLptr[j] = tmp;
                ce += cRLptr[j];
                        
                if(thinit == 1)
                    {
                    drawTZRlptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 54), j));
                    drawTZRlMptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 56), j));
                    cRLLptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 34), j));
                            
                    for(ii=0;ii<dPRlptr[j];++ii)
                        {
                        drawTZRlptr[ii+itrthin*dPRlptr[j]] = DRAWERS1Lptr[ii];
                        drawTZRlMptr[ii+itrthin*dPRlptr[j]] = DRAWERSM1Lptr[ii];
			if(all > 1)
                        	cRLLptr[ii+itrthin*dPRlptr[j]] = ce - nilor1Lptr[ii];
			else
				cRLLptr[ii+itrthin*dPRlptr[j]] = ce;
                        }
                    }
                }    
            }
        else
            {
            lambd = sig2/sigma2Rptr[j];
            s2tmp = 0.0;                    
            for(ii=0;ii<BRptr[j];++ii)
                {
                helpme = 1/(ZZptr[ii] + lambd*sptr[ii]);
                mur = helpme * nilorptr[ii];
                sigmaR = helpme * sig2;
                sigmaR = sqrt(sigmaR);
                GetRNGstate();
                drawtmp = rnorm(mur,sigmaR);
                PutRNGstate();
                DRAWERSptr[ii] = drawtmp;
                DRAWERSMptr[ii] = mur;
                s2tmp += (drawtmp*drawtmp)*sptr[ii];
                } 
                              
            // retransform parameters    
            cvp1(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 7), j), 2), VECTOR_ELT(VECTOR_ELT(ins, 59), j), VECTOR_ELT(VECTOR_ELT(ins, 10), j));        
            }
   
        // sampling random variance
        hypb = hypbt2 + s2tmp/2;
        sigma2Rptr[j] = 1/rand_gamma(haRnewptr[j],hypb);  
        
        ce -= cRptr[j];
        cRptr[j] = 0;
        for(ii=0;ii<BRptr[j];++ii)
            {
            cRptr[j] += cRMptr[ii]*drawBetaRptr[ii];
            }
        ce += cRptr[j];

        if(bcheckptr[j] > 0)
            {
            cvp1(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(ins, 7), j), 0), VECTOR_ELT(VECTOR_ELT(ins, 59), j), VECTOR_ELT(ins, 68));
                    
            mymean = 0.0;
            for(ii=0;ii<n;++ii)
                {
                mymean += fptr[ii];
                }
            mymean = mymean/n;
                    
            for(ii=0;ii<n;++ii)
                {
                fRptr[n*j + ii] = fptr[ii] - mymean;
                etaptr[ii] += fRptr[n*j + ii];
                }
            }
        else
            {
	    if(rccptr[j] > 0)
	    	{
		mymean = 0.0;
                for(ii=0;ii<BRptr[j];++ii)
                	{
                        mymean += drawBetaRptr[ii];
                        }
                mymean = mymean/BRptr[j];

                for(ii=0;ii<n;++ii)
                	{
                        // fptr[ii] = drawBetaRptr[rindptr[ii]];
                        // fRptr[n*j + ii] = fptr[ii];
			fRptr[n*j + ii] = drawBetaRptr[rindptr[ii]] - mymean;
                        etaptr[ii] += fRptr[n*j + ii];
                        }     
		} 
	   else
		{
                for(ii=0;ii<n;++ii)
                	{
			fRptr[n*j + ii] = drawBetaRptr[rindptr[ii]];
                        etaptr[ii] += fRptr[n*j + ii];
                        }  
		}                           
            }
                    
        if(thinit == 1)
            {                                        
            drawBptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 28), j));
            drawBMptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 58), j));
                                
            S2Rptr[j+itrthin*r] = sigma2Rptr[j];
	    if(all > 1)
            	cRMSptr[j+itrthin*r] = ce - cRptr[j];
	    else
		cRMSptr[j+itrthin*r] = ce;
            
            if(hspecptr[j] > 0)
                {
                stetaRptr = REAL(VECTOR_ELT(VECTOR_ELT(ins, 73), j));
                for(ii=0;ii<BRptr[j];++ii)
                    {
                    drawBptr[ii+itrthin*BRptr[j]] = DRAWERSptr[ii]; // - etaRptr[ii]; // drawBetaRptr[ii];
                    drawBMptr[ii+itrthin*BRptr[j]] = DRAWERSMptr[ii];
		    stetaRptr[ii+itrthin*BRptr[j]] = DRAWERSptr[ii] - etaRptr[ii];
                    }                          
                }
            else
                {
                for(ii=0;ii<BRptr[j];++ii)
                    {
                    drawBptr[ii+itrthin*BRptr[j]] = DRAWERSptr[ii]; // drawBetaRptr[ii];
                    drawBMptr[ii+itrthin*BRptr[j]] = DRAWERSMptr[ii];
                    }   
                }                    
            }
        }
    }
  

SEXP sampling(SEXP N, 
              SEXP ITER,
              SEXP thinsteps,
              SEXP response,
              SEXP sigma2, 
              SEXP K, 
              SEXP eta, 
              SEXP flin, 
              SEXP Ixx, 
              SEXP tildex, 
              SEXP ttildex, 
              SEXP hypbsigma2, 
              SEXP hypasigma2, 
              SEXP S2, 
              SEXP drawtb, 
              SEXP mdrawtb, 
              SEXP nlin,
              SEXP tbrq,
              SEXP bm,
              SEXP cL,
              SEXP CE,
              SEXP smstuff,
              SEXP mustuff,
              SEXP M,
              SEXP rastuff,
              SEXP R,
              SEXP PROB,
              SEXP low,
              SEXP high,
	      SEXP ITCPT,
	      SEXP ALL
              )
    {
    int i,ii,j,jj,d,dd,nProtected = 0,mcmc = 1000,thinit,itrthin = 0;
    
    PROTECT(N = coerceVector(N, INTSXP));
    ++nProtected;
    PROTECT(ITER = coerceVector(ITER, INTSXP));
    ++nProtected;
    PROTECT(K = coerceVector(K, INTSXP));
    ++nProtected;
    PROTECT(M = coerceVector(M, INTSXP));
    ++nProtected;
    PROTECT(nlin = coerceVector(nlin, INTSXP));
    ++nProtected;
    PROTECT(R = coerceVector(R, INTSXP));
    ++nProtected;
    
    const int iter = INTEGER(ITER)[0]; 
    const int k = INTEGER(K)[0];
    const int n = INTEGER(N)[0];
    const int nl = INTEGER(nlin)[0];
    const int m = INTEGER(M)[0];
    const int r = INTEGER(R)[0];
    const int itcpt = INTEGER(ITCPT)[0];
    const int all = INTEGER(ALL)[0];
    
    Rboolean prob = LOGICAL(PROB)[0];
    
    // setup smooth stuff
    double lambd,hypb,mus,drawtmp,hypbt2,tmp,s2tmp;
    double *fpsptr,*tau2ptr,*lambdaptr,*drawtbptr,*mdrawtbptr,*fsptr;
    double *niloptr,*drawssptr,*mussptr,*drawtgptr,*mdrawtgptr;
    double *hypatau2ptr,helpme,*cKptr,*ncoefptr,*csKptr,*sceptr,*spptr;
    int *dPptr,*cmeptr,*newsmptr,*lbdcheckptr;
    hypbt2 = 0.0001;
    hypb = 0;
    fpsptr =0;
    tau2ptr = 0;
    lambdaptr = 0;
    fsptr = 0;
    dPptr = 0;
    hypatau2ptr = 0;
    csKptr = 0;
    cmeptr = 0;
    sceptr = 0;
    lambd = 0;
    drawssptr = 0;
    mussptr = 0;
    drawtgptr = 0;
    mdrawtgptr = 0;
    tmp = 0;
    spptr = 0;
    newsmptr = 0;
    s2tmp = 0;
    
    SEXP fps,tildez,tau2,OldK,ttildez,hyperbtau2;
    SEXP hypatau2,lambda,dP,drawtg,mdrawtg;
    SEXP cme,cmcheck,cK,csK,sce;
    SEXP fs,nilo,drawss,muss,NCOEF,newsm,lbdcheck;
    fps = R_NilValue;
    tildez = R_NilValue;
    tau2 = R_NilValue;
    OldK = R_NilValue;
    ttildez = R_NilValue;
    hyperbtau2 = R_NilValue;
    hypatau2 = R_NilValue;
    lambda = R_NilValue;
    dP = R_NilValue;
    drawtg = R_NilValue;
    mdrawtg = R_NilValue;
    cme = R_NilValue;
    cmcheck = R_NilValue;
    cK = R_NilValue;
    csK = R_NilValue;
    sce = R_NilValue;
    fs = R_NilValue;
    nilo = R_NilValue;
    drawss = R_NilValue;
    muss = R_NilValue;
    NCOEF = R_NilValue;
    newsm = R_NilValue;
    lbdcheck = R_NilValue;
  
    if(k > 0)
        {
        PROTECT(fps = VECTOR_ELT(smstuff, 0));
        ++nProtected;  
        PROTECT(tildez = VECTOR_ELT(smstuff, 1));  
        ++nProtected;
        PROTECT(tau2 = VECTOR_ELT(smstuff, 2)); 
        ++nProtected;
        PROTECT(OldK = VECTOR_ELT(smstuff, 3)); 
        ++nProtected;
        PROTECT(ttildez = VECTOR_ELT(smstuff, 4)); 
        ++nProtected;
        PROTECT(hyperbtau2 = VECTOR_ELT(smstuff, 5)); 
        ++nProtected;
        PROTECT(hypatau2 = VECTOR_ELT(smstuff, 6));  
        ++nProtected;
        PROTECT(lambda = VECTOR_ELT(smstuff, 7));
        ++nProtected;
        PROTECT(dP = VECTOR_ELT(smstuff, 8)); 
        ++nProtected;
        PROTECT(drawtg = VECTOR_ELT(smstuff, 9)); 
        ++nProtected;
        PROTECT(mdrawtg = VECTOR_ELT(smstuff, 10));
        ++nProtected; 
        PROTECT(cme = VECTOR_ELT(smstuff, 11));
        ++nProtected;
        PROTECT(cmcheck = VECTOR_ELT(smstuff, 12));
        ++nProtected;
        PROTECT(sce = VECTOR_ELT(smstuff, 13));
        ++nProtected;
        PROTECT(cK = VECTOR_ELT(smstuff, 14));
        ++nProtected;
        PROTECT(csK = VECTOR_ELT(smstuff, 15));
        ++nProtected;
        PROTECT(newsm = VECTOR_ELT(smstuff, 16));
        ++nProtected;
	PROTECT(lbdcheck = VECTOR_ELT(smstuff, 17));
	++nProtected;
        
        PROTECT(fs = allocVector(VECSXP,k));
        ++nProtected;
    
        PROTECT(drawss = allocVector(VECSXP,k));
        ++nProtected;
    
        PROTECT(muss = allocVector(VECSXP,k));
        ++nProtected;
        
        PROTECT(nilo = allocVector(VECSXP,k));
        ++nProtected;
    
        PROTECT(NCOEF = allocVector(VECSXP,k));
        ++nProtected;
         
        for(j=0;j<k;++j)
            {       
            SEXP tmp1;
            PROTECT(tmp1 = allocVector(REALSXP,INTEGER(dP)[j]));
            SET_VECTOR_ELT(nilo, j, tmp1);
            UNPROTECT(1);
        
            SEXP tmp2;
            PROTECT(tmp2 = allocVector(REALSXP,INTEGER(dP)[j]));
            SET_VECTOR_ELT(drawss, j, tmp2);
            UNPROTECT(1);
        
            SEXP tmp3;
            PROTECT(tmp3 = allocVector(REALSXP,INTEGER(dP)[j]));
            SET_VECTOR_ELT(muss, j, tmp3);
            UNPROTECT(1);
        
            SEXP tmp4;
            PROTECT(tmp4 = allocVector(REALSXP,INTEGER(dP)[j]));
            SET_VECTOR_ELT(NCOEF, j, tmp4);
            UNPROTECT(1);
            
            SET_VECTOR_ELT(fs, j, allocVector(REALSXP, nrows(VECTOR_ELT(VECTOR_ELT(tildez, j), 0))));
            SET_VECTOR_ELT(newsm, j, coerceVector(VECTOR_ELT(newsm, j), INTSXP));
            }   
        hypbt2 = REAL(hyperbtau2)[0];
        fpsptr = REAL(fps);
        tau2ptr = REAL(tau2);
        lambdaptr = REAL(lambda);
        dPptr = INTEGER(dP);
        hypatau2ptr = REAL(hypatau2);
        csKptr = REAL(csK);
        cmeptr = INTEGER(cme);
        sceptr = REAL(sce);
	lbdcheckptr = INTEGER(lbdcheck);
        }
        
    // setup multilevel stuff
    SEXP fM,spec,hypident,tau2M,nlm,mlambdacheck,tildeZmu;
    SEXP ind,dPM,cM,cMM,fpsM,mcentering,hyperaTau2newM,lambdaM;
    SEXP drawTildeGammaM,drawTGMM,csM;
    SEXP eok,respok,nilook,indok,drawssM,mussM,ncoefM,ttildeZmu,mlcheck;
    fM = R_NilValue;
    spec = R_NilValue;
    hypident = R_NilValue;
    tau2M = R_NilValue;
    nlm = R_NilValue;
    mlambdacheck = R_NilValue;
    tildeZmu = R_NilValue;
    ind = R_NilValue;
    dPM = R_NilValue;
    cM = R_NilValue;
    cMM = R_NilValue;
    fpsM = R_NilValue;
    mcentering = R_NilValue;
    hyperaTau2newM = R_NilValue;
    lambdaM = R_NilValue;
    drawTildeGammaM = R_NilValue;
    drawTGMM = R_NilValue;
    csM = R_NilValue;
    eok = R_NilValue;
    respok = R_NilValue;
    nilook = R_NilValue;
    indok = R_NilValue;
    drawssM = R_NilValue;
    mussM = R_NilValue;
    ncoefM = R_NilValue;
    ttildeZmu = R_NilValue;
    mlcheck = R_NilValue;
    double *fMptr,*tau2Mptr,sum2 = 0,*mlambdacheckptr;
    double *eokptr,*respokptr,*nilookptr,*drawssMptr,*mussMptr,*ncoefMptr;
    double *cMptr,*cMcKptr,*fpsMptr,*hyperaTau2newMptr;
    double *drawTildeGammaMptr,*drawTGMMptr,*lambdaMptr,*csMptr,*sMptr;
    int *specptr,*nlmptr,*dPMptr,*indokptr,*indptr,*mlcheckptr;
    fMptr = 0;
    specptr = 0;
    nlmptr = 0;
    mlambdacheckptr = 0;
    dPMptr = 0;
    indokptr = 0;
    eokptr = 0;
    respokptr = 0;
    indptr = 0;
    nilookptr = 0;
    drawssMptr = 0;
    mussMptr = 0;
    ncoefMptr = 0;
    cMptr = 0;
    cMcKptr = 0;
    fpsMptr = 0;
    hyperaTau2newMptr = 0;
    drawTildeGammaMptr = 0;
    drawTGMMptr = 0;
    lambdaMptr = 0;
    csMptr = 0;
    mlcheckptr = 0;
    sMptr = 0;
   
    if(m > 0)
        {
        PROTECT(fM = VECTOR_ELT(mustuff, 0));
        ++nProtected;
        PROTECT(spec = VECTOR_ELT(mustuff, 1));
        ++nProtected;
        PROTECT(hypident = VECTOR_ELT(mustuff, 2));
        ++nProtected;
        PROTECT(tau2M = VECTOR_ELT(mustuff, 3));
        ++nProtected;
        PROTECT(nlm = VECTOR_ELT(mustuff, 4));
        ++nProtected;
        PROTECT(mlambdacheck = VECTOR_ELT(mustuff, 5));
        ++nProtected;
        PROTECT(tildeZmu = VECTOR_ELT(mustuff, 6));
        ++nProtected;
        PROTECT(ind = VECTOR_ELT(mustuff, 7));
        ++nProtected;
        PROTECT(dPM = VECTOR_ELT(mustuff, 8));
        ++nProtected;
        PROTECT(cM = VECTOR_ELT(mustuff, 9));
        ++nProtected;
        PROTECT(cMM = VECTOR_ELT(mustuff, 10));
        ++nProtected;
        PROTECT(fpsM = VECTOR_ELT(mustuff, 11));
        ++nProtected;
        PROTECT(mcentering = VECTOR_ELT(mustuff, 12));
        ++nProtected;
        PROTECT(hyperaTau2newM = VECTOR_ELT(mustuff, 13));
        ++nProtected;
        PROTECT(lambdaM = VECTOR_ELT(mustuff, 14));
        ++nProtected;
        PROTECT(drawTildeGammaM = VECTOR_ELT(mustuff, 15));
        ++nProtected;
        PROTECT(drawTGMM = VECTOR_ELT(mustuff, 16));
        ++nProtected;
        PROTECT(csM = VECTOR_ELT(mustuff, 17));
        ++nProtected;
        PROTECT(ttildeZmu = VECTOR_ELT(mustuff, 18));
        ++nProtected;
        PROTECT(mlcheck = VECTOR_ELT(mustuff, 19));
        ++nProtected;

        PROTECT(fM = coerceVector(fM, REALSXP));
        ++nProtected;
        PROTECT(spec = coerceVector(spec, INTSXP));
        ++nProtected;   
        PROTECT(nlm = coerceVector(nlm, INTSXP));
        ++nProtected; 
        fMptr = REAL(fM); 
        specptr = INTEGER(spec);
        nlmptr = INTEGER(nlm);
        
        PROTECT(eok = allocVector(VECSXP, m));
        ++nProtected; 
        PROTECT(respok = allocVector(VECSXP, m));
        ++nProtected;
        PROTECT(nilook = allocVector(VECSXP, m));
        ++nProtected;
        PROTECT(indok = allocVector(VECSXP, m));
        ++nProtected;
        PROTECT(drawssM = allocVector(VECSXP, m));
        ++nProtected;
        PROTECT(mussM = allocVector(VECSXP, m));
        ++nProtected;
        PROTECT(ncoefM = allocVector(VECSXP, m));
        ++nProtected;
                
        for(j=0;j<m;++j)
            {   
            SEXP tmp6;
            PROTECT(tmp6 = allocVector(VECSXP, nlmptr[j]));   
            SEXP tmp7;
            PROTECT(tmp7 = allocVector(VECSXP, nlmptr[j]));      
            SEXP tmp8;
            PROTECT(tmp8 = allocVector(VECSXP, nlmptr[j]));
            SEXP tmp12;
            PROTECT(tmp12 = allocVector(INTSXP, nlmptr[j]));
            
            SET_VECTOR_ELT(dPM, j, coerceVector(VECTOR_ELT(dPM, j), INTSXP));
            
            SET_VECTOR_ELT(drawssM, j, allocVector(REALSXP, INTEGER(VECTOR_ELT(dPM, j))[0]));
            SET_VECTOR_ELT(mussM, j, allocVector(REALSXP, INTEGER(VECTOR_ELT(dPM, j))[0]));
            SET_VECTOR_ELT(ncoefM, j, allocVector(REALSXP, INTEGER(VECTOR_ELT(dPM, j))[0]));
           
            for(jj=0;jj<nlmptr[j];++jj)
                {
                int lobsm = length(VECTOR_ELT(VECTOR_ELT(ind, j), jj));
                INTEGER(tmp12)[jj] = lobsm;
                
                SEXP tmp9;
                PROTECT(tmp9 = allocVector(REALSXP,lobsm));
                SET_VECTOR_ELT(tmp6, jj, tmp9);
                UNPROTECT(1);
                
                SEXP tmp10;
                PROTECT(tmp10 = allocVector(REALSXP,lobsm));
                SET_VECTOR_ELT(tmp7, jj, tmp10);
                UNPROTECT(1);

                SEXP tmp11;
                PROTECT(tmp11 = allocVector(REALSXP,INTEGER(VECTOR_ELT(dPM, j))[jj]));
                SET_VECTOR_ELT(tmp8, jj, tmp11);
                UNPROTECT(1);
                
                SET_VECTOR_ELT(VECTOR_ELT(ind, j), jj, coerceVector(VECTOR_ELT(VECTOR_ELT(ind, j), jj), INTSXP));

                for(ii=0;ii<lobsm;++ii)
                    INTEGER(VECTOR_ELT(VECTOR_ELT(ind, j), jj))[ii] -= 1;
                }
           
            SET_VECTOR_ELT(nilook, j, tmp8);
            UNPROTECT(1);        
            SET_VECTOR_ELT(respok, j, tmp7);
            UNPROTECT(1);       
            SET_VECTOR_ELT(eok, j, tmp6);
            UNPROTECT(1); 
            SET_VECTOR_ELT(indok, j, tmp12);
            UNPROTECT(1); 
            }
        }
           
    // setup random stuff
    SEXP fR,hspec,rwcheck,sigma2R,Z,etaR,BR,drawBetaR,cR,cRM,RRi,ins;
    SEXP metaR,RRs,fr,la,tau2R,tildeZra,OTzRa,dPR,ttildeZra,center;
    SEXP hyperaTau2newR,cRS,cRMSS,lambdaR,drawTildeGammaR,drawTildeGammaRM;
    SEXP cRSS,newR,RRl,frl,tZRl,ttZRl,cRL,cRML,drawTZRl,drawTZRlM,cRLL;
    SEXP haRnew,drawB,drawBM,S2R,cRMS;
    SEXP DRAWERS,DRAWERSM,nilor;
    SEXP DRAWERS1,DRAWERSM1,nilor1;
    SEXP DRAWERS1L,DRAWERSM1L,nilor1L,dPRl,fRtmp;
    SEXP nilore,rind,ZZ,bcheck,stetaR,rcc;

    fR = R_NilValue;
    hspec = R_NilValue;
    rwcheck = R_NilValue;
    sigma2R = R_NilValue;
    Z = R_NilValue;
    etaR = R_NilValue;
    BR = R_NilValue;
    drawBetaR = R_NilValue;
    cR = R_NilValue;
    cRM = R_NilValue;
    RRi = R_NilValue;
    ins = R_NilValue;
    metaR = R_NilValue;
    RRs = R_NilValue;
    fr = R_NilValue;
    la = R_NilValue;
    tau2R = R_NilValue;
    tildeZra = R_NilValue;
    OTzRa = R_NilValue;
    dPR = R_NilValue;
    ttildeZra = R_NilValue;
    center = R_NilValue;
    hyperaTau2newR = R_NilValue;
    cRS = R_NilValue;
    cRMSS = R_NilValue;
    lambdaR = R_NilValue;
    drawTildeGammaR = R_NilValue;
    drawTildeGammaRM = R_NilValue;
    cRSS = R_NilValue;
    newR = R_NilValue;
    RRl = R_NilValue;
    frl = R_NilValue;
    tZRl = R_NilValue;
    ttZRl = R_NilValue;
    cRL = R_NilValue;
    cRML = R_NilValue;
    drawTZRl = R_NilValue;
    drawTZRlM = R_NilValue;
    cRLL = R_NilValue;
    haRnew = R_NilValue;
    drawB = R_NilValue;
    drawBM = R_NilValue;
    S2R = R_NilValue;
    cRMS = R_NilValue;
    DRAWERS = R_NilValue;
    DRAWERSM = R_NilValue;
    nilor = R_NilValue;
    DRAWERS1 = R_NilValue;
    DRAWERSM1 = R_NilValue;
    nilor1 = R_NilValue;
    DRAWERS1L = R_NilValue;
    DRAWERSM1L = R_NilValue;
    nilor1L = R_NilValue;
    dPRl = R_NilValue;
    fRtmp = R_NilValue;
    nilore = R_NilValue;
    rind = R_NilValue;
    ZZ = R_NilValue;
    bcheck = R_NilValue;
    stetaR = R_NilValue;
    rcc = R_NilValue;
    
    double *fRptr,*sigma2Rptr,*sptr,*etaRptr;
    double *DRAWERSptr,*DRAWERSMptr,*haRnewptr,*drawBptr,*drawBMptr;
    double mur = 0,*S2Rptr,*nilorptr,*drawBetaRptr,*cRptr,*cRMptr,*cRMSptr;
    double *frptr,*tau2Rptr,*hyperaTau2newRptr;
    double *DRAWERS1ptr,*DRAWERSM1ptr,*nilor1ptr,*s1ptr;
    double *drawTildeGammaRptr,*drawTildeGammaRMptr,*frlptr;
    double *metaRptr,*cRLptr,*cRMLptr,*drawTZRlptr,*drawTZRlMptr,*cRLLptr;
    double *DRAWERS1Lptr,*DRAWERSM1Lptr,*nilor1Lptr,sigmaR,*cRSptr,*cRMSSptr,*fRtmpptr;
    double *ZZptr,*stetaRptr;
    fRptr = 0;
    sigma2Rptr = 0;
    sptr = 0;
    etaRptr = 0;
    DRAWERSptr = 0;
    DRAWERSMptr = 0;
    haRnewptr = 0;
    drawBptr = 0;
    drawBMptr = 0;
    S2Rptr = 0;
    nilorptr = 0;
    drawBetaRptr = 0;
    cRptr = 0;
    cRMptr = 0;
    cRMSptr = 0;
    frptr = 0;
    tau2Rptr = 0;
    DRAWERS1ptr = 0;
    DRAWERSM1ptr = 0;
    nilor1ptr = 0;
    s1ptr = 0;
    hyperaTau2newRptr = 0;
    drawTildeGammaRptr = 0;
    drawTildeGammaRMptr = 0;
    frlptr = 0;
    metaRptr = 0;
    cRLptr = 0;
    cRMLptr = 0;
    drawTZRlptr = 0;
    drawTZRlMptr = 0;
    cRLLptr = 0;
    DRAWERS1Lptr = 0;
    DRAWERSM1Lptr = 0;
    nilor1Lptr = 0;
    cRSptr = 0;
    cRMSSptr = 0;
    fRtmpptr = 0;
    ZZptr = 0;
    stetaRptr = 0;

    int *hspecptr,*BRptr,*RRsptr,*RRlptr,*dPRptr,*centerptr,*dPRlptr;    
    int *RRiptr,*newRptr,*rindptr,*bcheckptr,*rccptr;
    hspecptr = 0;
    BRptr = 0;
    RRsptr = 0;
    RRlptr = 0;
    dPRptr = 0;
    centerptr = 0;
    dPRlptr = 0;
    RRiptr = 0;
    newRptr = 0;
    rindptr = 0;
    bcheckptr = 0;
    rccptr = 0;
        
    if(r > 0)
        {
        PROTECT(fR = VECTOR_ELT(rastuff, 0));
        ++nProtected;
        PROTECT(hspec = coerceVector(VECTOR_ELT(rastuff, 1), INTSXP));
        ++nProtected;
        PROTECT(rwcheck = VECTOR_ELT(rastuff, 2));
        ++nProtected;
        PROTECT(sigma2R = VECTOR_ELT(rastuff, 3));
        ++nProtected;
        PROTECT(Z = VECTOR_ELT(rastuff, 4));
        ++nProtected;
        PROTECT(etaR = VECTOR_ELT(rastuff, 5));
        ++nProtected;
        PROTECT(BR = coerceVector(VECTOR_ELT(rastuff, 6), INTSXP));
        ++nProtected;
        PROTECT(drawBetaR = VECTOR_ELT(rastuff, 7));
        ++nProtected;
        PROTECT(cR = VECTOR_ELT(rastuff, 8));
        ++nProtected;
        PROTECT(cRM = VECTOR_ELT(rastuff, 9));
        ++nProtected;
        PROTECT(RRi = coerceVector(VECTOR_ELT(rastuff, 10), INTSXP));
        ++nProtected;
        PROTECT(ins = VECTOR_ELT(rastuff, 11));
        ++nProtected;
        PROTECT(metaR = VECTOR_ELT(rastuff, 12));
        ++nProtected;
        PROTECT(RRs = coerceVector(VECTOR_ELT(rastuff, 13), INTSXP));
        ++nProtected;
        PROTECT(fr = VECTOR_ELT(rastuff, 14));
        ++nProtected;
        PROTECT(la = VECTOR_ELT(rastuff, 15));
        ++nProtected;
        PROTECT(tau2R = VECTOR_ELT(rastuff, 16));
        ++nProtected;
        PROTECT(tildeZra = VECTOR_ELT(rastuff, 17));
        ++nProtected;
        PROTECT(OTzRa = VECTOR_ELT(rastuff, 18));
        ++nProtected;
        PROTECT(dPR = VECTOR_ELT(rastuff, 19));
        ++nProtected;
        PROTECT(ttildeZra = VECTOR_ELT(rastuff, 20));
        ++nProtected;
        PROTECT(center = VECTOR_ELT(rastuff, 21));
        ++nProtected;
        PROTECT(hyperaTau2newR = VECTOR_ELT(rastuff, 22));
        ++nProtected;
        PROTECT(cRS = VECTOR_ELT(rastuff, 23));
        ++nProtected;
        PROTECT(cRMSS = VECTOR_ELT(rastuff, 24));
        ++nProtected;
        PROTECT(lambdaR = VECTOR_ELT(rastuff, 25));
        ++nProtected;
        PROTECT(drawTildeGammaR = VECTOR_ELT(rastuff, 26));
        ++nProtected;
        PROTECT(drawTildeGammaRM = VECTOR_ELT(rastuff, 27));
        ++nProtected;
        PROTECT(cRSS = VECTOR_ELT(rastuff, 28));
        ++nProtected;
        PROTECT(newR = VECTOR_ELT(rastuff, 29));
        ++nProtected;
        PROTECT(RRl = coerceVector(VECTOR_ELT(rastuff, 30), INTSXP));
        ++nProtected;
        PROTECT(frl = VECTOR_ELT(rastuff, 31));
        ++nProtected;
        PROTECT(tZRl = VECTOR_ELT(rastuff, 32));
        ++nProtected;
        PROTECT(ttZRl = VECTOR_ELT(rastuff, 33));
        ++nProtected;
        PROTECT(cRL = VECTOR_ELT(rastuff, 34));
        ++nProtected;
        PROTECT(cRML = VECTOR_ELT(rastuff, 35));
        ++nProtected;
        PROTECT(drawTZRl = VECTOR_ELT(rastuff, 36));
        ++nProtected;
        PROTECT(drawTZRlM = VECTOR_ELT(rastuff, 37));
        ++nProtected;
        PROTECT(cRLL = VECTOR_ELT(rastuff, 38));
        ++nProtected;
        PROTECT(haRnew = VECTOR_ELT(rastuff, 39));
        ++nProtected;
        PROTECT(drawB = VECTOR_ELT(rastuff, 40));
        ++nProtected;
        PROTECT(drawBM = VECTOR_ELT(rastuff, 41));
        ++nProtected;
        PROTECT(S2R = VECTOR_ELT(rastuff, 42));
        ++nProtected;
        PROTECT(cRMS = VECTOR_ELT(rastuff, 43));        
        ++nProtected;
        PROTECT(dPRl = coerceVector(VECTOR_ELT(rastuff, 44), INTSXP));
        ++nProtected;
        PROTECT(fRtmp = VECTOR_ELT(rastuff, 45));
        ++nProtected;
        PROTECT(rind = VECTOR_ELT(rastuff, 46));
        ++nProtected;
        PROTECT(ZZ = VECTOR_ELT(rastuff, 47));
        ++nProtected;
        PROTECT(bcheck = coerceVector(VECTOR_ELT(rastuff, 48), INTSXP));
        ++nProtected;
	PROTECT(stetaR = VECTOR_ELT(rastuff, 49));
        ++nProtected;
	PROTECT(rcc = VECTOR_ELT(rastuff, 50));
        ++nProtected;
        
        fRptr = REAL(fR);
        sigma2Rptr = REAL(sigma2R);
        S2Rptr = REAL(S2R);
        cRMSptr = REAL(cRMS);
        cRLptr = REAL(cRL);
        haRnewptr = REAL(haRnew);
        cRptr = REAL(cR);              
        hspecptr = INTEGER(hspec);
        BRptr = INTEGER(BR);
        RRsptr = INTEGER(RRs);
        RRlptr = INTEGER(RRl);
        dPRlptr = INTEGER(dPRl);
        RRiptr = INTEGER(RRi);
        bcheckptr = INTEGER(bcheck);
	rccptr = INTEGER(rcc);
                
        PROTECT(DRAWERS = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(DRAWERSM = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(nilor = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(nilore = allocVector(VECSXP, r));
        ++nProtected;
        
        PROTECT(DRAWERS1 = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(DRAWERSM1 = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(nilor1 = allocVector(VECSXP, r));
        ++nProtected;
        
        PROTECT(DRAWERS1L = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(DRAWERSM1L = allocVector(VECSXP, r));
        ++nProtected;
        PROTECT(nilor1L = allocVector(VECSXP, r));
        ++nProtected;
        
        for(j=0;j<r;++j)
            {
            SET_VECTOR_ELT(DRAWERS, j, allocVector(REALSXP,BRptr[j]));
            SET_VECTOR_ELT(DRAWERSM, j, allocVector(REALSXP,BRptr[j]));
            SET_VECTOR_ELT(nilor, j, allocVector(REALSXP,BRptr[j]));
            SET_VECTOR_ELT(nilore, j, allocVector(REALSXP,BRptr[j]));
            
            if(RRsptr[j] > 0)
                {
                SET_VECTOR_ELT(DRAWERS1, j, allocVector(VECSXP, RRsptr[j]));
                SET_VECTOR_ELT(DRAWERSM1, j, allocVector(VECSXP, RRsptr[j]));
                SET_VECTOR_ELT(nilor1, j, allocVector(VECSXP, RRsptr[j]));
                SET_VECTOR_ELT(dPR, j, coerceVector(VECTOR_ELT(dPR, j), INTSXP));
                SET_VECTOR_ELT(center, j, coerceVector(VECTOR_ELT(center, j), INTSXP));
            
                for(jj=0;jj<RRsptr[j];++jj)
                    {
                    SET_VECTOR_ELT(VECTOR_ELT(DRAWERS1, j), jj, allocVector(REALSXP,INTEGER(VECTOR_ELT(dPR, j))[jj]));
                    SET_VECTOR_ELT(VECTOR_ELT(DRAWERSM1, j), jj, allocVector(REALSXP,INTEGER(VECTOR_ELT(dPR, j))[jj]));
                    SET_VECTOR_ELT(VECTOR_ELT(nilor1, j), jj, allocVector(REALSXP,INTEGER(VECTOR_ELT(dPR, j))[jj]));
                    }
                }
            if(RRlptr[j] > 0)
                {
                SET_VECTOR_ELT(DRAWERS1L, j, allocVector(REALSXP,RRlptr[j]));
                SET_VECTOR_ELT(DRAWERSM1L, j, allocVector(REALSXP,RRlptr[j]));
                SET_VECTOR_ELT(nilor1L, j, allocVector(REALSXP,RRlptr[j]));
                }
            }
        }
             
    // setup the rest
    SEXP e,mu,draws,tx,meta,bdraws,f,origresponse;
    double *eptr,*metaptr,mymean,*fptr;
   
    PROTECT(meta = allocVector(REALSXP,n));        
    ++nProtected; 
    metaptr = REAL(meta);
    PROTECT(e = allocVector(REALSXP,n));
    ++nProtected;
    PROTECT(f = allocVector(REALSXP,n));
    ++nProtected;
    PROTECT(origresponse = allocVector(REALSXP,n));        
    ++nProtected; 
    eptr = REAL(e);
    fptr = REAL(f);
    for(i=0;i<n;++i)
        {
        eptr[i] = 0.0;
        metaptr[i] = 0.0;
	REAL(origresponse)[i] = REAL(response)[i];
        }  
    
    PROTECT(mu = allocVector(REALSXP,nl));
    ++nProtected;
    
    PROTECT(draws = allocVector(REALSXP,nl));
    ++nProtected;
    
    PROTECT(tx = VECTOR_ELT(tildex, 0));
    ++nProtected;
    
    PROTECT(bdraws = allocVector(REALSXP,nl));
    ++nProtected;
        
    double sumtmp=0.0,hyp1,hyp2,diceta=0.0,sigma=0.001;
    hyp1 = REAL(hypasigma2)[0];
    hyp2 = REAL(hypbsigma2)[0];
    
    for(i=0;i<length(thinsteps);++i)
        {
        INTEGER(thinsteps)[i] = INTEGER(thinsteps)[i] - 1;
        }
        
    time_t start,end;
    double dif,timeout;    
    double *sigma2ptr,*S2ptr,sig2,sig,*etaptr,*flinptr,*responseptr;
    double ce = 0.0,bs = 0.0;
    double *bdrawsptr,*cLptr,*CEptr,*muptr,*drawsptr;
    sigma2ptr = REAL(sigma2);
    S2ptr = REAL(S2);
    sig2 = REAL(sigma2)[0];
    sig = sqrt(sig2);
    etaptr = REAL(eta);
    flinptr = REAL(flin);
    responseptr = REAL(response);  
    bdrawsptr = REAL(bdraws);
    drawsptr = REAL(draws);
    drawtbptr = REAL(drawtb);
    mdrawtbptr = REAL(mdrawtb);
    muptr = REAL(mu);
    cLptr = REAL(cL);
    CEptr = REAL(CE);   
    start = time(NULL);
    
    int mcount = 0, maxit = 20;
    double mcsmu = 0;
        
    Rprintf("Starting\n");    
     
    GetRNGstate();
    for(i=0;i<iter;++i)
        {
        thinit = thincheck(thinsteps,i);  
        
        if(prob)
            {
            vtnorm(low,high,eta,1,response,n);
            }
        
        // sampling for smooth effects
        if(k > 0)
            { 
            for(j=0;j<k;++j)
                {     
                for(ii=0;ii<n;++ii)
                    {
                    etaptr[ii] -= fpsptr[n*j + ii];
                    eptr[ii] = responseptr[ii] - etaptr[ii];
                    }

                niloptr = REAL(VECTOR_ELT(nilo, j));
                drawssptr = REAL(VECTOR_ELT(drawss, j));
                mussptr = REAL(VECTOR_ELT(muss, j));
                cKptr = REAL(VECTOR_ELT(cK, j));
                ncoefptr = REAL(VECTOR_ELT(NCOEF, j));
                spptr = REAL(VECTOR_ELT(VECTOR_ELT(tildez,j), 1));
                fsptr = REAL(VECTOR_ELT(fs, j));
                newsmptr = INTEGER(VECTOR_ELT(newsm , j));
        
		if(lbdcheckptr[j]<1) 
                	lambd = sig2/tau2ptr[j];
		else
			lambd = REAL(VECTOR_ELT(VECTOR_ELT(tildez,j), 4))[0];               
		if(INTEGER(VECTOR_ELT(VECTOR_ELT(tildez,j), 3))[0] > 0)
			{
			monysamp(maxit, 
              			 VECTOR_ELT(VECTOR_ELT(tildez,j), 4), 
              			 VECTOR_ELT(VECTOR_ELT(tildez,j), 1), 
              			 lambd, 
              			 dPptr[j],
              			 VECTOR_ELT(VECTOR_ELT(tildez,j), 0), 
              			 e, 
              			 VECTOR_ELT(muss, j), 
              			 VECTOR_ELT(drawss, j), 
              		  	 sumtmp, 
              			 sig2, 
              			 INTEGER(VECTOR_ELT(VECTOR_ELT(tildez,j), 3))[0]);

/*			sigmat(VECTOR_ELT(VECTOR_ELT(tildez,j), 4),VECTOR_ELT(VECTOR_ELT(tildez,j), 1),lambd,VECTOR_ELT(sigmats,j),dPptr[j]);*/
/*			cvp2(VECTOR_ELT(VECTOR_ELT(tildez,j), 0),e,VECTOR_ELT(nilo, j));*/
/*			sigmatsptr = REAL(VECTOR_ELT(sigmats,j));*/

/*			mcount = 0;*/
/*			while(mcount < maxit)*/
/*				{*/
/*				for(ii=0;ii<dPptr[j];++ii)*/
/*					{*/
/*					mcsmu = 0;*/

/*					if(ii>0)*/
/*						for(dd=0;dd<ii;++dd)*/
/*							mcsmu += drawssptr[dd]*sigmatsptr[ii + dPptr[j]*dd];*/
/*					if(ii<(dPptr[j] - 1))*/
/*						for(dd=(ii+1);dd<dPptr[j];++dd)*/
/*							mcsmu += drawssptr[dd]*sigmatsptr[ii + dPptr[j]*dd];*/
/*				*/
/*					mcsmu = (niloptr[ii] - mcsmu)/sigmatsptr[ii + dPptr[j]*ii];*/

/*					if(INTEGER(VECTOR_ELT(VECTOR_ELT(tildez,j), 3))[0] == 1)*/
/*						{*/
/*						if(ii==0)*/
/*							drawssptr[ii] = trunc_normal2(-1000,drawssptr[1],mcsmu,sqrt(sig2/sigmatsptr[ii + dPptr[j]*ii]));*/
/*						if(ii==(dPptr[j]-1))*/
/*							drawssptr[ii] = trunc_normal2(drawssptr[ii-1],1000,mcsmu,sqrt(sig2/sigmatsptr[ii + dPptr[j]*ii]));*/
/*						if(ii>0 && ii < (dPptr[j]-1)) */
/*							drawssptr[ii] = trunc_normal2(drawssptr[ii-1],drawssptr[ii+1],mcsmu,sqrt(sig2/sigmatsptr[ii + dPptr[j]*ii]));*/
/*						}*/
/*					else*/
/*						{*/
/*						if(ii==0)*/
/*							drawssptr[ii] = trunc_normal2(drawssptr[1],1000,mcsmu,sqrt(sig2/sigmatsptr[ii + dPptr[j]*ii]));*/
/*						if(ii==(dPptr[j]-1))*/
/*							drawssptr[ii] = trunc_normal2(-1000,drawssptr[ii-1],mcsmu,sqrt(sig2/sigmatsptr[ii + dPptr[j]*ii]));*/
/*						if(ii>0 && ii < (dPptr[j]-1)) */
/*							drawssptr[ii] = trunc_normal2(drawssptr[ii+1],drawssptr[ii-1],mcsmu,sqrt(sig2/sigmatsptr[ii + dPptr[j]*ii]));*/
/*						}*/
/*					}*/
/*				mcount++;*/
/*				}*/

/*			sumtmp = sums(VECTOR_ELT(drawss, j), VECTOR_ELT(VECTOR_ELT(tildez,j), 1));*/
			}
		else
			{                               
                	cvp2(VECTOR_ELT(OldK, j),e,VECTOR_ELT(nilo, j));               
                	sumtmp = 0.0;
                	for(ii=0;ii<dPptr[j];++ii)
                    		{
                    		helpme = 1/(1+lambd*spptr[ii]);    
                    		mus = helpme * niloptr[ii];
                    		sigma = helpme * sig2;
                    		sigma = sqrt(sigma);
				GetRNGstate();
                    		drawtmp = rnorm(mus, sigma);
				PutRNGstate();
                    		drawssptr[ii] = drawtmp;
                    		mussptr[ii] = mus;
                    		sumtmp += (drawtmp*drawtmp)*spptr[ii];
                    		} 
			}
// check !!!!!!!!!
/*		if(cmeptr[j] == 1)*/
/*			rtr(spptr,drawssptr,lambd,dPptr[j],sumtmp);*/

                hypb = hypbt2 + sumtmp/2;
                tau2ptr[j] = 1/rand_gamma(hypatau2ptr[j],hypb);


                ce -= sceptr[j];
                sceptr[j] = 0;
                if(LOGICAL(VECTOR_ELT(VECTOR_ELT(cmcheck, j), 0))[0])
                    {
                    cvp2(VECTOR_ELT(VECTOR_ELT(cmcheck, j), 1),VECTOR_ELT(drawss, j),VECTOR_ELT(NCOEF, j));
                    for(ii=0;ii<dPptr[j];++ii)
                        sceptr[j] += cKptr[ii]*ncoefptr[ii];
                    }
                else
                    {
                    for(ii=0;ii<dPptr[j];++ii)
                        sceptr[j] += cKptr[ii]*drawssptr[ii];
                    }
                ce += sceptr[j];
        
                cvp2(VECTOR_ELT(ttildez, j),VECTOR_ELT(drawss, j),VECTOR_ELT(fs, j));
                
                if(cmeptr[j] == 1)
                    {
                    gcenter(VECTOR_ELT(fs, j)); 
                    for(ii=0;ii<n;++ii)
                        {
                        fpsptr[n*j + ii] = fsptr[newsmptr[ii]];
                        etaptr[ii] += fpsptr[n*j + ii];
                        }
                    }
                else
                    {
                    for(ii=0;ii<n;++ii)
                        {
                        fpsptr[n*j + ii] = fsptr[newsmptr[ii]];
                        etaptr[ii] += fsptr[ii];
                        }                    
                    }
/*                for(ii=0;ii<n;++ii)*/
/*                	{*/
/*                        fpsptr[n*j + ii] = fsptr[newsmptr[ii]];*/
/*                        etaptr[ii] += fsptr[ii];*/
/*                        }   */
                
                if(thinit == 1)
                    {
                    drawtgptr = REAL(VECTOR_ELT(drawtg, j));
                    mdrawtgptr = REAL(VECTOR_ELT(mdrawtg, j));
                    
                    lambdaptr[j+itrthin*k] = tau2ptr[j];
		    if(all > 1)
                    	csKptr[j+itrthin*k] = ce - sceptr[j];
		    else
                    	csKptr[j+itrthin*k] = sceptr[j];

                    for(ii=0;ii<dPptr[j];++ii)
                        {
                        drawtgptr[ii+itrthin*dPptr[j]] = drawssptr[ii];
                        mdrawtgptr[ii+itrthin*dPptr[j]] = mussptr[ii];
                        }
                    }
                }
            }
            
        // sampling multilevel effects
        if(m > 0)
            {
            for(j=0;j<m;++j)
                {                   
                if(specptr[j] == 2)
                    {
                    tau2Mptr = REAL(VECTOR_ELT(tau2M, j));
                    mlambdacheckptr = REAL(VECTOR_ELT(mlambdacheck, j));
                    dPMptr = INTEGER(VECTOR_ELT(dPM, j));
                    cMptr = REAL(VECTOR_ELT(cM, j));
                    drawssMptr = REAL(VECTOR_ELT(drawssM, j));
                    mussMptr = REAL(VECTOR_ELT(mussM, j));
                    ncoefMptr = REAL(VECTOR_ELT(ncoefM, j));
                    indokptr = INTEGER(VECTOR_ELT(indok, j));
                    hyperaTau2newMptr = REAL(VECTOR_ELT(hyperaTau2newM, j));
                    lambdaMptr = REAL(VECTOR_ELT(lambdaM, j));
                    csMptr = REAL(VECTOR_ELT(csM, j));
                    mlcheckptr = INTEGER(VECTOR_ELT(mlcheck, j));
                    
                    if(LOGICAL(hypident)[j])
                        {
                        sum2 = 0;
                        lambd = sig2/tau2Mptr[0];
                        }
                        
                    for(jj=0;jj<nlmptr[j];++jj)
                        {
                        if(!LOGICAL(hypident)[j])
                            {
                            if(mlcheckptr[jj] < 1)
                                {
                                lambd = sig2/tau2Mptr[jj];
                                }
                            else
                                {
                                lambd = mlambdacheckptr[jj];
                                }
                            }
                            
                        eokptr = REAL(VECTOR_ELT(VECTOR_ELT(eok, j), jj));   
                        respokptr = REAL(VECTOR_ELT(VECTOR_ELT(respok, j), jj));
                        indptr = INTEGER(VECTOR_ELT(VECTOR_ELT(ind, j), jj));
                        nilookptr = REAL(VECTOR_ELT(VECTOR_ELT(nilook, j), jj));
               
                        for(ii=0;ii<indokptr[jj];++ii)
                            {
                            etaptr[indptr[ii]] -= fMptr[n*j + indptr[ii]];
                            eokptr[ii] = responseptr[indptr[ii]] - etaptr[indptr[ii]];
                            }
                        
                        if(INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 3))[0] > 0)
				{
				monysamp(maxit, 
              			 	 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 4), 
              			 	 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 1), 
              			 	 lambd, 
              			 	 dPMptr[jj],
              			 	 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 0), 
              			 	 VECTOR_ELT(VECTOR_ELT(eok, j), jj), 
              			 	 VECTOR_ELT(mussM, j), 
              			 	 VECTOR_ELT(drawssM, j), 
              		  	 	 sumtmp, 
              			 	 sig2, 
              			 	 INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 3))[0]);
				}
			else
				{
                        	cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 0), 
				     VECTOR_ELT(VECTOR_ELT(eok, j), jj), VECTOR_ELT(VECTOR_ELT(nilook, j), jj));  
                        
                        	sMptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 1));
                        	sumtmp = 0.0;
                        	for(ii=0;ii<dPMptr[jj];++ii)
                             		{
                             		helpme = 1/(1+lambd*sMptr[ii]);    
                             		mus = helpme * nilookptr[ii];
                             		sigma = helpme * sig2;
                             		sigma = sqrt(sigma);
			     		GetRNGstate();
                             		drawtmp = rnorm(mus, sigma);
			     		PutRNGstate();
                             		drawssMptr[ii] = drawtmp;
                             		mussMptr[ii] = mus;
                             		sumtmp += (drawtmp*drawtmp)*sMptr[ii];
                             		} 
				}
  
                        ce -= cMptr[jj];
                        cMptr[jj] = 0;
                        cMcKptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cMM, j), jj), 1));
                    
                        if(LOGICAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cMM, j), jj), 0))[0])
                            {
                            cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cMM, j), jj), 2), VECTOR_ELT(drawssM, j), VECTOR_ELT(ncoefM, j));
                            for(ii=0;ii<dPMptr[jj];++ii)
                                cMptr[jj] += cMcKptr[ii]*ncoefMptr[ii];
                            }
                        else
                            {
                            for(ii=0;ii<dPMptr[jj];++ii)
                                cMptr[jj] += cMcKptr[ii]*drawssMptr[ii];
                            }
                        ce += cMptr[jj];    

/*			if(LOGICAL(VECTOR_ELT(mcentering, j))[jj])*/
/*				rtr(sMptr,drawssMptr,lambd,dPMptr[jj],sumtmp);*/
                        
                        cvp2(VECTOR_ELT(VECTOR_ELT(ttildeZmu, j), jj), VECTOR_ELT(drawssM, j), VECTOR_ELT(VECTOR_ELT(fpsM, j), jj));
                        
                        fpsMptr = REAL(VECTOR_ELT(VECTOR_ELT(fpsM, j), jj));
                        if(LOGICAL(VECTOR_ELT(mcentering, j))[jj])
                            {
                            //gcenter(VECTOR_ELT(VECTOR_ELT(fpsM, j), jj));
                            mymean = 0.0;
                            for(ii=0;ii<indokptr[jj];++ii)
                                {
                                mymean += fpsMptr[ii];
                                }
                            mymean = mymean/indokptr[jj];
                            for(ii=0;ii<indokptr[jj];++ii)
                                {
                                fMptr[n*j + indptr[ii]] = fpsMptr[ii] - mymean;
                                etaptr[indptr[ii]] += fMptr[n*j + indptr[ii]];
                                }
                            }
                        else
                            {
                            for(ii=0;ii<indokptr[jj];++ii)
                                {
                                fMptr[n*j + indptr[ii]] = fpsMptr[ii];
                                etaptr[indptr[ii]] += fMptr[n*j + indptr[ii]];
                                }                            
                            }

/*                        for(ii=0;ii<indokptr[jj];++ii)*/
/*                        	{*/
/*                                fMptr[n*j + indptr[ii]] = fpsMptr[ii];*/
/*                                etaptr[indptr[ii]] += fMptr[n*j + indptr[ii]];*/
/*                                }   */
                        if(!LOGICAL(hypident)[j])
                            {    
                            hypb = hypbt2 + sumtmp/2;
                            tau2Mptr[jj] = 1/rand_gamma(hyperaTau2newMptr[jj],hypb); 
                            }
                        else
                            {
                            sum2 += sumtmp;
                            }
                      
                        if(thinit == 1)
                            {
                            drawTildeGammaMptr = REAL(VECTOR_ELT(VECTOR_ELT(drawTildeGammaM, j), jj));
                            drawTGMMptr = REAL(VECTOR_ELT(VECTOR_ELT(drawTGMM, j), jj));
                    
                            if(LOGICAL(hypident)[j])
                                {
                                lambdaMptr[jj+itrthin*nlmptr[j]] = tau2Mptr[0];
                                }
                            else
                                {
                                lambdaMptr[jj+itrthin*nlmptr[j]] = tau2Mptr[jj];
                                }
                            
			    if(all > 1)    
                            	csMptr[jj+itrthin*nlmptr[j]] = ce - cMptr[jj];
			    else
                    		csMptr[jj+itrthin*nlmptr[j]] = ce;
                            for(ii=0;ii<dPMptr[jj];++ii)
                                {
                                drawTildeGammaMptr[ii+itrthin*dPMptr[jj]] = drawssMptr[ii];
                                drawTGMMptr[ii+itrthin*dPMptr[jj]] = mussMptr[ii];
                                }
                            }          
                        }
                    if(LOGICAL(hypident)[j])
                        {
                        hypb = hypbt2 + sum2/2;
                        tau2Mptr[0] = 1/rand_gamma(hyperaTau2newMptr[0],hypb);
                        }
                    }
                if(specptr[j] == 1)
                    {
                    cMptr = REAL(VECTOR_ELT(cM, j));
                    drawssMptr = REAL(VECTOR_ELT(drawssM, j));
                    mussMptr = REAL(VECTOR_ELT(mussM, j));
                    ncoefMptr = REAL(VECTOR_ELT(ncoefM, j));
                    indokptr = INTEGER(VECTOR_ELT(indok, j));
                    dPMptr = INTEGER(VECTOR_ELT(dPM, j));

                    for(jj=0;jj<nlmptr[j];++jj)
                        {                   
                        eokptr = REAL(VECTOR_ELT(VECTOR_ELT(eok, j), jj));   
                        respokptr = REAL(VECTOR_ELT(VECTOR_ELT(respok, j), jj));
                        indptr = INTEGER(VECTOR_ELT(VECTOR_ELT(ind, j), jj));
                        nilookptr = REAL(VECTOR_ELT(VECTOR_ELT(nilook, j), jj));
                        csMptr = REAL(VECTOR_ELT(VECTOR_ELT(csM, j), jj));
                         
                        for(ii=0;ii<indokptr[jj];++ii)
                           {
                           etaptr[indptr[ii]] -= fMptr[n*j + indptr[ii]];
                           eokptr[ii] = responseptr[indptr[ii]] - etaptr[indptr[ii]];
                           }

                        cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZmu, j), jj), 0), VECTOR_ELT(VECTOR_ELT(eok, j), jj), VECTOR_ELT(VECTOR_ELT(nilook, j), jj));  
                        for(ii=0;ii<dPMptr[jj];++ii)
                            {
                            drawtmp = rnorm(nilookptr[ii], sig);
                            drawssMptr[ii] = drawtmp;
                            mussMptr[ii] = nilookptr[ii];
                            } 

                        ce -= cMptr[jj];
                        cMptr[jj] = 0;
                        cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cMM, j), jj), 0), VECTOR_ELT(drawssM, j), VECTOR_ELT(VECTOR_ELT(nilook, j), jj));
                        for(ii=0;ii<dPMptr[jj];++ii)
                            cMptr[jj] += nilookptr[ii];
                        ce += cMptr[jj];

                        cvp2(VECTOR_ELT(VECTOR_ELT(ttildeZmu, j), jj), VECTOR_ELT(drawssM, j), VECTOR_ELT(VECTOR_ELT(eok, j), jj));
                        for(ii=0;ii<indokptr[jj];++ii)
                            {
                            fMptr[n*j + indptr[ii]] = eokptr[ii];
                            etaptr[indptr[ii]] += fMptr[n*j + indptr[ii]];
                            }

                        if(thinit == 1)
                            {
                            drawTildeGammaMptr = REAL(VECTOR_ELT(VECTOR_ELT(drawTildeGammaM, j), jj));
                            drawTGMMptr = REAL(VECTOR_ELT(VECTOR_ELT(drawTGMM, j), jj));
                            
                            for(ii=0;ii<dPMptr[jj];++ii)
                                {
                                drawTildeGammaMptr[ii+itrthin*dPMptr[jj]] = drawssMptr[ii];
                                drawTGMMptr[ii+itrthin*dPMptr[jj]] = mussMptr[ii];
				if(all > 1)
                                	csMptr[ii+itrthin*dPMptr[jj]] = ce - nilookptr[ii];
				else
					csMptr[ii+itrthin*dPMptr[jj]] = ce;
                                }
                            } 
                        }
                    }
                }
            }
            
        // sampling random effects
        if(r > 0)
            {
            for(j=0;j<r;++j)
                {  
                sptr = REAL(VECTOR_ELT(VECTOR_ELT(Z, j), 1));
                DRAWERSptr = REAL(VECTOR_ELT(DRAWERS, j));
                DRAWERSMptr = REAL(VECTOR_ELT(DRAWERSM, j));
                drawBetaRptr = REAL(VECTOR_ELT(drawBetaR, j));
                nilorptr = REAL(VECTOR_ELT(nilor, j));
                cRMptr = REAL(VECTOR_ELT(cRM, j));
                rindptr = INTEGER(VECTOR_ELT(rind, j));
                ZZptr = REAL(VECTOR_ELT(ZZ, j));
 
                for(ii=0;ii<n;++ii)
                    {
                    etaptr[ii] -= fRptr[n*j + ii];
                    eptr[ii] = responseptr[ii] - etaptr[ii];
                    }

                cvp2(VECTOR_ELT(VECTOR_ELT(Z, j), 0), e, VECTOR_ELT(nilor, j));
                                
                if(hspecptr[j] > 0)
                    {
                    etaRptr = REAL(VECTOR_ELT(etaR, j));
                    metaRptr = REAL(VECTOR_ELT(metaR, j));

		    s2tmp = 0.0;
                    
                    if(LOGICAL(rwcheck)[j])
                        {
                        lambd = sig2/sigma2Rptr[j];
                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            helpme = 1/(ZZptr[ii] + lambd*sptr[ii]);
                            mur = helpme*(nilorptr[ii] + lambd*etaRptr[ii]);
                            sigmaR = helpme * sig2;
                            sigmaR = sqrt(sigmaR);
			    GetRNGstate();
                            drawtmp = rnorm(mur,sigmaR);
			    PutRNGstate();
			    tmp = drawtmp - etaRptr[ii];
			    s2tmp += tmp*tmp;
                            DRAWERSptr[ii] = drawtmp;
                            DRAWERSMptr[ii] = mur;
                            }
                        }
                    else
                        {
                        lambd = sig2/sigma2Rptr[j];
                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            tmp = lambd*sptr[ii];
                            helpme = 1/(ZZptr[ii]  + tmp);
                            mur = helpme*(nilorptr[ii] + tmp*etaRptr[ii]);
                            sigmaR = helpme * sig2;
                            sigmaR = sqrt(sigmaR);
	                    GetRNGstate();
                            drawtmp = rnorm(mur,sigmaR);
                            PutRNGstate();
			    tmp = drawtmp - etaRptr[ii];
			    s2tmp += tmp*tmp;
                            DRAWERSptr[ii] = drawtmp;
                            DRAWERSMptr[ii] = mur;
                            }
                        }

                    
/*                    for(ii=0;ii<BRptr[j];++ii)*/
/*                        {*/
/*                        tmp = drawBetaRptr[ii] - etaRptr[ii];*/
/*                        s2tmp += tmp*tmp;*/
/*                        }  */
                        
                    // retransform parameters   
                    cvp1(VECTOR_ELT(VECTOR_ELT(Z, j), 2), VECTOR_ELT(DRAWERS, j), VECTOR_ELT(drawBetaR, j)); 

                    if(RRiptr[j] > 0)
                        {
                        rsampling(VECTOR_ELT(ins, j),VECTOR_ELT(drawBetaR, j),VECTOR_ELT(etaR, j),VECTOR_ELT(nilore, j),
				  BRptr[j],sigma2Rptr[j],hypbt2,hypb,ce,thinit,itrthin,maxit,all); 
                        }
                    if(RRsptr[j] > 0)
                        {
                        tau2Rptr = REAL(VECTOR_ELT(tau2R, j));
                        dPRptr = INTEGER(VECTOR_ELT(dPR, j));
                        hyperaTau2newRptr = REAL(VECTOR_ELT(hyperaTau2newR, j));
                        centerptr = INTEGER(VECTOR_ELT(center, j));
                        cRSptr = REAL(VECTOR_ELT(cRS, j));
                        
                        for(jj=0;jj<RRsptr[j];++jj)
                            {
                            s1ptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZra, j), jj), 1));
                            frptr = REAL(VECTOR_ELT(VECTOR_ELT(fr, j), jj));
                            DRAWERS1ptr = REAL(VECTOR_ELT(VECTOR_ELT(DRAWERS1, j), jj));
                            DRAWERSM1ptr = REAL(VECTOR_ELT(VECTOR_ELT(DRAWERSM1, j), jj));
                            nilor1ptr = REAL(VECTOR_ELT(VECTOR_ELT(nilor1, j), jj));   
                            newRptr = INTEGER(VECTOR_ELT(VECTOR_ELT(newR, j), jj));
                            fRtmpptr = REAL(VECTOR_ELT(VECTOR_ELT(fRtmp, j), jj));
                            
                            for(ii=0;ii<BRptr[j];++ii)
                                {
                                etaRptr[ii] -= frptr[ii];
                                nilorptr[ii] = drawBetaRptr[ii] - etaRptr[ii];
                                }
                                
                            if(LOGICAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(la, j), jj), 0))[0])
                                {
                                lambd = sigma2Rptr[j]/tau2Rptr[jj];
                                }
                            else
                                {
                                lambd = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(la, j), jj), 1))[0];
                                }

			    if(INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZra, j), jj), 3))[0] > 0)
			    	{
				monysamp(maxit, 
              			 	 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZra, j), jj), 4), 
              			 	 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZra, j), jj), 1), 
              			 	 lambd, 
              			 	 dPRptr[jj],
              			 	 VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZra, j), jj), 0), 
              			 	 VECTOR_ELT(VECTOR_ELT(nilor1, j), jj), 
              			 	 VECTOR_ELT(VECTOR_ELT(DRAWERSM1, j), jj), 
              			 	 VECTOR_ELT(VECTOR_ELT(DRAWERS1, j), jj), 
              		  	 	 sumtmp, 
              			 	 sigmaR, 
              			 	 INTEGER(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(tildeZra, j), jj), 3))[0]);
				}
			    else
				{           
                            	cvp2(VECTOR_ELT(VECTOR_ELT(OTzRa, j), jj), VECTOR_ELT(nilor, j), VECTOR_ELT(VECTOR_ELT(nilor1, j), jj));
                                
                            	sumtmp = 0.0;                    
                            	for(ii=0;ii<dPRptr[jj];++ii)
                                	{
                                	helpme = 1/(1 + lambd*s1ptr[ii]);
                                	mur = helpme * nilor1ptr[ii];
                                	sigmaR = helpme * sigma2Rptr[j];
                                	sigmaR = sqrt(sigmaR);
	                        	GetRNGstate();
                                	drawtmp = rnorm(mur,sigmaR);
					PutRNGstate();
                                	DRAWERS1ptr[ii] = drawtmp;
                                	DRAWERSM1ptr[ii] = mur;
                                	sumtmp += (drawtmp*drawtmp)*s1ptr[ii];
                                	}
				}

/*			    if(centerptr[jj] > 0)*/
/*			    	rtr(s1ptr,DRAWERS1ptr,lambd,dPRptr[jj],sumtmp);*/

                            // sampling smooth variance
                            hypb = hypbt2 + sumtmp/2;
                            tau2Rptr[jj] = 1/rand_gamma(hyperaTau2newRptr[jj],hypb);
                           
                            cvp2(VECTOR_ELT(VECTOR_ELT(ttildeZra, j), jj), VECTOR_ELT(VECTOR_ELT(DRAWERS1, j), jj), VECTOR_ELT(VECTOR_ELT(fRtmp, j), jj));
 
                            if(centerptr[jj] > 0)
                                {
                                gcenter(VECTOR_ELT(VECTOR_ELT(fRtmp, j), jj));
                                for(ii=0;ii<BRptr[j];++ii)
                                   {
                                   frptr[ii] = fRtmpptr[newRptr[ii]];
                                   etaRptr[ii] += frptr[ii];
                                   metaRptr[ii] += etaRptr[ii];
                                   }
                                }
                            else
                                {
                                for(ii=0;ii<BRptr[j];++ii)
                                    {
                                    frptr[ii] = fRtmpptr[newRptr[ii]];
                                    etaRptr[ii] += frptr[ii];
                                    metaRptr[ii] += etaRptr[ii];
                                    }                                
                                }

/*                            for(ii=0;ii<BRptr[j];++ii)*/
/*                            	{*/
/*                                frptr[ii] = fRtmpptr[newRptr[ii]];*/
/*                                etaRptr[ii] += frptr[ii];*/
/*                                metaRptr[ii] += etaRptr[ii];*/
/*                                }     */
                   
                            ce -= cRSptr[jj];
                            cRSptr[jj] = 0;
                            cRMSSptr = REAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cRMSS, j), jj), 0));

                            if(LOGICAL(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cRMSS, j), jj), 1))[0])
                                {
                                cvp2(VECTOR_ELT(VECTOR_ELT(VECTOR_ELT(cRMSS, j), jj), 2),
				     VECTOR_ELT(VECTOR_ELT(DRAWERS1, j), jj),VECTOR_ELT(VECTOR_ELT(nilor1, j), jj));

                                for(ii=0;ii<dPRptr[jj];++ii)
                                    cRSptr[jj] += cRMSSptr[ii]*nilor1ptr[ii];
                                }
                            else
                                {
                                for(ii=0;ii<dPRptr[jj];++ii)
                                    cRSptr[jj] += cRMSSptr[ii]*DRAWERS1ptr[ii];
                                }
                            ce += cRSptr[jj];
                             
                            if(thinit == 1)
                                {
                                REAL(VECTOR_ELT(VECTOR_ELT(lambdaR, j), jj))[itrthin] = tau2Rptr[jj];
				if(all > 1)
                                	REAL(VECTOR_ELT(VECTOR_ELT(cRSS, j), jj))[itrthin] = ce - cRSptr[jj];
				else
                                	REAL(VECTOR_ELT(VECTOR_ELT(cRSS, j), jj))[itrthin] = ce;
                                drawTildeGammaRptr = REAL(VECTOR_ELT(VECTOR_ELT(drawTildeGammaR, j), jj));
                                drawTildeGammaRMptr = REAL(VECTOR_ELT(VECTOR_ELT(drawTildeGammaRM, j), jj));
                                
                                for(ii=0;ii<dPRptr[jj];++ii)
                                    {
                                    drawTildeGammaRptr[ii+itrthin*dPRptr[jj]] = DRAWERS1ptr[ii];
                                    drawTildeGammaRMptr[ii+itrthin*dPRptr[jj]] = DRAWERSM1ptr[ii];
                                    }
                                }
                            }
                        }
                    if(RRlptr[j] > 0)
                        {
                        frlptr = REAL(VECTOR_ELT(frl, j));
                        DRAWERS1Lptr = REAL(VECTOR_ELT(DRAWERS1L, j));
                        DRAWERSM1Lptr = REAL(VECTOR_ELT(DRAWERSM1L, j));
                        cRMLptr = REAL(VECTOR_ELT(VECTOR_ELT(cRML, j), 1));
                        nilor1Lptr = REAL(VECTOR_ELT(nilor1L, j));
                        
                        for(ii=0;ii<BRptr[j];++ii)
                                {
                                etaRptr[ii] -= frlptr[ii];
                                nilorptr[ii] = drawBetaRptr[ii] - etaRptr[ii];
                                }
                        cvp2(VECTOR_ELT(tZRl, j), VECTOR_ELT(nilor, j), VECTOR_ELT(DRAWERSM1L, j));
                        
                        sigmaR = sqrt(sigma2Rptr[j]);
                        for(ii=0;ii<dPRlptr[j];++ii)
                            {
			    GetRNGstate();
                            DRAWERS1Lptr[ii] = rnorm(DRAWERSM1Lptr[ii], sigmaR);
			    PutRNGstate();
                            }
                        cvp1(VECTOR_ELT(tZRl, j),VECTOR_ELT(DRAWERS1L, j),VECTOR_ELT(frl, j));

                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            etaRptr[ii] += frlptr[ii];
                            metaRptr[ii] += etaRptr[ii];
                            }

                        ce -= cRLptr[j];
                        cvp1(VECTOR_ELT(VECTOR_ELT(cRML, j), 0),VECTOR_ELT(DRAWERS1L, j),VECTOR_ELT(nilor1L, j));
                        tmp = 0;
                        for(ii=0;ii<dPRlptr[j];++ii)
                            {
                            nilor1Lptr[ii] *= cRMLptr[ii];
                            tmp += nilor1Lptr[ii];
                            }
                        cRLptr[j] = tmp;
                        ce += cRLptr[j];
                        
                        if(thinit == 1)
                            {
                            drawTZRlptr = REAL(VECTOR_ELT(drawTZRl, j));
                            drawTZRlMptr = REAL(VECTOR_ELT(drawTZRlM, j));
                            cRLLptr = REAL(VECTOR_ELT(cRLL, j));
                            
                            for(ii=0;ii<dPRlptr[j];++ii)
                                {
                                drawTZRlptr[ii+itrthin*dPRlptr[j]] = DRAWERS1Lptr[ii];
                                drawTZRlMptr[ii+itrthin*dPRlptr[j]] = DRAWERSM1Lptr[ii];
				if(all > 1)
                                	cRLLptr[ii+itrthin*dPRlptr[j]] = ce - nilor1Lptr[ii];
				else
					cRLLptr[ii+itrthin*dPRlptr[j]] = ce;
                                }
                            }
                        }   
                    }
                else
                    {       
                    lambd = sig2/sigma2Rptr[j];
                    s2tmp = 0.0;             
                    for(ii=0;ii<BRptr[j];++ii)
                        {
                        helpme = 1/(ZZptr[ii]  + lambd*sptr[ii]);
                        mur = helpme * nilorptr[ii];
                        sigmaR = helpme * sig2;
                        sigmaR = sqrt(sigmaR);
                        GetRNGstate();
                        drawtmp = rnorm(mur,sigmaR);
                        PutRNGstate();
                        DRAWERSptr[ii] = drawtmp;
                        DRAWERSMptr[ii] = mur;
                        s2tmp += (drawtmp*drawtmp)*sptr[ii];
                        } 
                              
                    // retransform parameters    
                    cvp1(VECTOR_ELT(VECTOR_ELT(Z, j), 2), VECTOR_ELT(DRAWERS, j), VECTOR_ELT(drawBetaR, j));        
                    }  
                               
                // draws for random variance
                hypb = hypbt2 + s2tmp/2;
                sigma2Rptr[j] = 1/rand_gamma(haRnewptr[j],hypb);                   
                
                ce -= cRptr[j];
                cRptr[j] = 0;
                for(ii=0;ii<BRptr[j];++ii)
                    cRptr[j] += cRMptr[ii]*drawBetaRptr[ii];
                ce += cRptr[j];

                if(bcheckptr[j] > 0)
                    {
                    cvp1(VECTOR_ELT(VECTOR_ELT(Z, j), 0), VECTOR_ELT(DRAWERS, j), f);
                    
                    mymean = 0.0;
                    for(ii=0;ii<n;++ii)
                        {
                        mymean += fptr[ii];
                        }
                    mymean = mymean/n;
                    
                    for(ii=0;ii<n;++ii)
                        {
                        fRptr[n*j + ii] = fptr[ii] - mymean;
                        etaptr[ii] += fRptr[n*j + ii];
                        }
                    }
                else
                    {
		    if(rccptr[j] > 0)
			{
		    	mymean = 0.0;
                    	for(ii=0;ii<BRptr[j];++ii)
                        	{
                        	mymean += drawBetaRptr[ii];
                        	}
                    	mymean = mymean/BRptr[j];

                    	for(ii=0;ii<n;++ii)
                        	{
                        	// fptr[ii] = drawBetaRptr[rindptr[ii]];
                        	// fRptr[n*j + ii] = fptr[ii];
				fRptr[n*j + ii] = drawBetaRptr[rindptr[ii]] - mymean;
                        	etaptr[ii] += fRptr[n*j + ii];
                        	}     
			} 
		    else
			{
                    	for(ii=0;ii<n;++ii)
                        	{
				fRptr[n*j + ii] = drawBetaRptr[rindptr[ii]];
                        	etaptr[ii] += fRptr[n*j + ii];
                        	}  
			}              
                    }
                    
                if(thinit == 1)
                    {                                        
                    drawBptr = REAL(VECTOR_ELT(drawB, j));
                    drawBMptr = REAL(VECTOR_ELT(drawBM, j));
                    
                    S2Rptr[j+itrthin*r] = sigma2Rptr[j];
		    if(all > 1)
                    	cRMSptr[j+itrthin*r] = ce - cRptr[j];
		    else
                        cRMSptr[j+itrthin*r] = ce;
                    if(hspecptr[j] > 0)
                        {
			stetaRptr = REAL(VECTOR_ELT(stetaR, j));
                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            drawBptr[ii+itrthin*BRptr[j]] = DRAWERSptr[ii]; // - etaRptr[ii]; // drawBetaRptr[ii];
                            drawBMptr[ii+itrthin*BRptr[j]] = DRAWERSMptr[ii];
			    stetaRptr[ii+itrthin*BRptr[j]] = DRAWERSptr[ii] - etaRptr[ii];
                            }                          
                        }
                    else
                        {
                        for(ii=0;ii<BRptr[j];++ii)
                            {
                            drawBptr[ii+itrthin*BRptr[j]] = DRAWERSptr[ii]; // drawBetaRptr[ii];
                            drawBMptr[ii+itrthin*BRptr[j]] = DRAWERSMptr[ii];
                            }   
                        }                 
                    }
                }
            }

        // sampling linear effects
	if(itcpt > 0 || nl > 1)
		{
        	for(ii=0;ii<n;++ii)
            		{
            		etaptr[ii] = etaptr[ii] - flinptr[ii];
            		eptr[ii] = responseptr[ii] - etaptr[ii];
            		}
       
        	cvp2(tx,e,mu);
        	sigma = sqrt(sig2);
        	for(ii=0;ii<nl;++ii)
	    		{
	    		GetRNGstate();
            		drawsptr[ii] = rnorm(muptr[ii], sigma);
	    		PutRNGstate();
	    		}
        	cvp2(ttildex,draws,flin);
        
        	ce -= bs;
        	bs = 0;
        	cvp2(tbrq,draws,bdraws);
        	for(ii=0;ii<nl;++ii)
            		{
            		bs += bdrawsptr[ii];
           		}
        	ce += bs;
		}
        
        sumtmp = 0.0;
        for(ii=0;ii<n;++ii)
            {
            etaptr[ii] = etaptr[ii] + flinptr[ii];
            eptr[ii] = responseptr[ii] - etaptr[ii];
            sumtmp += eptr[ii]*eptr[ii];
            }
        
        // sampling the variance
        if(!prob)
            {
            hypb = hyp2 + sumtmp/2;        
            sig2 = 1/rand_gamma(hyp1,hypb);
            sig = sqrt(sig2);
            }
        
        // storing      
        if(thinit == 1)
            {
            S2ptr[itrthin] = sig2;
            for(ii=0;ii<nl;++ii)
                {
                drawtbptr[ii+itrthin*nl] = drawsptr[ii];
                mdrawtbptr[ii+itrthin*nl] = muptr[ii];
                if(ii > 0)
                    {
		    if(all > 1)
                    	cLptr[(ii-1)+itrthin*(nl-1)] = ce - bdrawsptr[ii];
		    else
			cLptr[(ii-1)+itrthin*(nl-1)] = ce;
                    }
                CEptr[itrthin] = ce;
                }
           
            for(ii=0;ii<n;++ii)
                {
                metaptr[ii] += etaptr[ii];
                }
            
            diceta += dev(response,eta,sig2,n);
             
            itrthin += 1;
            }
        if(i == mcmc)
            {
            Rprintf("MCMC iteration: %d\n",mcmc);
            mcmc += 1000;
            }
        if(i == 999)
            {
            end = time(NULL);
            dif = difftime(end,start);
            timeout = dif*iter/999;
            Rprintf("Approximate sampling time is: %g",timeout);
            Rprintf(" sec\n");
            
            if(timeout > 200)
                {
                Rprintf("Coffee break : ) \n\n"); 
                Rprintf("                 .Z:\n");                    
                Rprintf("                8MM.\n");                    
                Rprintf("              ,MMM2    ..\n");               
                Rprintf("              @MMM,    6MC\n");              
                Rprintf("           v  .MMMM    :MM.\n");             
                Rprintf("          ,MMc  $MMM    MM7\n");             
                Rprintf("           WMM$   MMM  .MM.\n");             
                Rprintf("            MMMc  @MM. MM0\n");              
                Rprintf("             MM# $MM; MM#\n");               
                Rprintf("              MM MM1  M#\n");                
                Rprintf("              n. oM\n");                     
                Rprintf("                   \n");                     
                Rprintf("                ..\n");                      
                Rprintf("         .AMMMMMMMMM;\n");                   
                Rprintf("      1MMMM@Ai.      cMMMMMMMM@:\n");        
                Rprintf("     @MMM.           MMM@0bMMMME\n");        
                Rprintf("     MMMM;                 MMMc vMMMMME\n"); 
                Rprintf("     .MMMMMM07i. . . .iZMMMM. 1MMM71MMME\n");
                Rprintf("       QMMMMMMMMMMMMMMMM@c  .MMM,   .MMM\n");
                Rprintf("        AMMM. .i:;ii..     nMMM     @MMZ\n");
                Rprintf("         SMMM             tMMM  @MMMMMM\n"); 
                Rprintf("          vMMM           ,MMM  YMMMMMM.\n"); 
                Rprintf("           9MMM          MMM     .;i\n");    
                Rprintf("     .;0MMMMMMM@        9MM  ,;;\n");        
                Rprintf("  c#MMMMMWt. .MMM      QMM  MMMMMMMi\n");    
                Rprintf("vMMMMMz        MMMMMMMMM@   ,0QEWMMM0\n");   
                Rprintf("MMMMMM1          .ctInv          MMMb\n");   
                Rprintf(" ;MMMMMMW7                    .BMMM7\n");    
                Rprintf("   .A@MMMMMMMMB9ti,.....i76@MMMM@i\n");      
                Rprintf("         :79WMMMMMMMMMMMMM@07.\n\n");  
                }
            }
        }
    PutRNGstate(); 

    if(prob)
	{
    	for(i=0;i<n;++i)
		REAL(response)[i] = REAL(origresponse)[i];
	}
    
    Rprintf("MCMC iteration: %d\n",iter);
   
    SEXP names;
    PROTECT(names = allocVector(STRSXP, 12));
    ++nProtected;
      
    SEXP out;  
    PROTECT(out = allocVector(VECSXP, 12));
    ++nProtected;
    
    SEXP dic;  
    PROTECT(dic = allocVector(REALSXP, 1));
    REAL(dic)[0] = diceta;
    ++nProtected;
    
    SEXP thin;  
    PROTECT(thin = allocVector(REALSXP, 1));
    REAL(thin)[0] = itrthin;
    ++nProtected;
    
    // return smooth, mu stuff and random stuff
    SEXP outsm,outmu,outra;
    outsm = R_NilValue;
    outmu = R_NilValue;
    outra = R_NilValue;
    
    SEXP namessm;
    PROTECT(namessm = allocVector(STRSXP, 4));
    ++nProtected;
    
    SEXP namesmu;
    PROTECT(namesmu = allocVector(STRSXP, 4));
    ++nProtected;
    
    SEXP namesra;
    PROTECT(namesra = allocVector(STRSXP, 14));
    ++nProtected;
    
    if(k > 0)
        {
        PROTECT(outsm = allocVector(VECSXP, 4));
        ++nProtected;
        
        SET_VECTOR_ELT(outsm, 0, lambda);
        SET_VECTOR_ELT(outsm, 1, drawtg);
        SET_VECTOR_ELT(outsm, 2, mdrawtg);
        SET_VECTOR_ELT(outsm, 3, csK);
        
        SET_STRING_ELT(namessm, 0, mkChar("lambda"));
        SET_STRING_ELT(namessm, 1, mkChar("drawtg"));
        SET_STRING_ELT(namessm, 2, mkChar("mdrawtg"));
        SET_STRING_ELT(namessm, 3, mkChar("csK"));
        
        setAttrib(outsm, R_NamesSymbol, namessm);        
        }
    if(m > 0)
        {
        PROTECT(outmu = allocVector(VECSXP, 4));
        ++nProtected;
        
        SET_VECTOR_ELT(outmu, 0, lambdaM);
        SET_VECTOR_ELT(outmu, 1, drawTildeGammaM);
        SET_VECTOR_ELT(outmu, 2, drawTGMM);
        SET_VECTOR_ELT(outmu, 3, csM);
        
        SET_STRING_ELT(namesmu, 0, mkChar("lambdaM"));
        SET_STRING_ELT(namesmu, 1, mkChar("drawTildeGammaM"));
        SET_STRING_ELT(namesmu, 2, mkChar("drawTGMM"));
        SET_STRING_ELT(namesmu, 3, mkChar("csM"));
        
        setAttrib(outmu, R_NamesSymbol, namesmu);
        }
    if(r > 0)
        {
        PROTECT(outra = allocVector(VECSXP, 14));
        ++nProtected;
        
        SET_VECTOR_ELT(outra, 0, drawB);
        SET_VECTOR_ELT(outra, 1, drawBM);
        SET_VECTOR_ELT(outra, 2, S2R);
        SET_VECTOR_ELT(outra, 3, cRMS);
        SET_VECTOR_ELT(outra, 4, lambdaR);
        SET_VECTOR_ELT(outra, 5, drawTildeGammaR);
        SET_VECTOR_ELT(outra, 6, drawTildeGammaRM);
        SET_VECTOR_ELT(outra, 7, cRSS);
        SET_VECTOR_ELT(outra, 8, drawTZRl);
        SET_VECTOR_ELT(outra, 9, drawTZRlM);
        SET_VECTOR_ELT(outra, 10, cRLL);
        SET_VECTOR_ELT(outra, 11, metaR);
        SET_VECTOR_ELT(outra, 12, ins);
        SET_VECTOR_ELT(outra, 13, stetaR);
        
        SET_STRING_ELT(namesra, 0, mkChar("drawB"));
        SET_STRING_ELT(namesra, 1, mkChar("drawBM"));
        SET_STRING_ELT(namesra, 2, mkChar("S2R"));
        SET_STRING_ELT(namesra, 3, mkChar("cRMS"));
        SET_STRING_ELT(namesra, 4, mkChar("lambdaR"));
        SET_STRING_ELT(namesra, 5, mkChar("drawTildeGammaR"));
        SET_STRING_ELT(namesra, 6, mkChar("drawTildeGammaRM"));
        SET_STRING_ELT(namesra, 7, mkChar("cRSS"));
        SET_STRING_ELT(namesra, 8, mkChar("drawTZRl"));
        SET_STRING_ELT(namesra, 9, mkChar("drawTZRlM"));
        SET_STRING_ELT(namesra, 10, mkChar("cRLL"));
        SET_STRING_ELT(namesra, 11, mkChar("metaR"));
        SET_STRING_ELT(namesra, 12, mkChar("ins"));
        SET_STRING_ELT(namesra, 13, mkChar("stetaR"));
        
        setAttrib(outra, R_NamesSymbol, namesra);
        }
        
    // return the rest
    SET_VECTOR_ELT(out, 0, S2);
    SET_VECTOR_ELT(out, 1, drawtb);
    SET_VECTOR_ELT(out, 2, mdrawtb);
    SET_VECTOR_ELT(out, 3, draws);
    SET_VECTOR_ELT(out, 4, meta);
    SET_VECTOR_ELT(out, 5, dic);
    SET_VECTOR_ELT(out, 6, thin);
    SET_VECTOR_ELT(out, 7, cL);
    SET_VECTOR_ELT(out, 8, CE);
    SET_VECTOR_ELT(out, 9, outsm);
    SET_VECTOR_ELT(out, 10, outmu);
    SET_VECTOR_ELT(out, 11, outra);
    
    SET_STRING_ELT(names, 0, mkChar("S2"));
    SET_STRING_ELT(names, 1, mkChar("lcoef"));
    SET_STRING_ELT(names, 2, mkChar("mlcoef"));
    SET_STRING_ELT(names, 3, mkChar("ldraws"));
    SET_STRING_ELT(names, 4, mkChar("meta"));
    SET_STRING_ELT(names, 5, mkChar("diceta"));
    SET_STRING_ELT(names, 6, mkChar("itrthin"));
    SET_STRING_ELT(names, 7, mkChar("cL"));
    SET_STRING_ELT(names, 8, mkChar("CE"));
    SET_STRING_ELT(names, 9, mkChar("outsm"));
    SET_STRING_ELT(names, 10, mkChar("outmu"));
    SET_STRING_ELT(names, 11, mkChar("outra"));
            
    setAttrib(out, R_NamesSymbol, names);
    
    UNPROTECT(nProtected);
    return out;   
    }
