#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>


SEXP tr(SEXP mat, SEXP umat, SEXP ind, SEXP ind2, SEXP uind, SEXP N, SEXP COL, SEXP NL)
    {
    int i,j,k,n,nl,col,one,two;
    
    n = INTEGER(N)[0];
    nl = INTEGER(NL)[0];
    col = INTEGER(COL)[0];
    
    int *iptr,*uiptr;
    iptr = INTEGER(ind);
    uiptr = INTEGER(uind);
    
    double *mptr,*umptr,*i2ptr;
    mptr = REAL(mat);
    umptr = REAL(umat);
    i2ptr = REAL(ind2);
    
    double s;
    
    for(j=0;j<nl;++j)
        {
        s = 0;
        two = uiptr[j];
        for(i=0;i<n;++i)
            {
            one = iptr[i];
            if(one == two)
                {
                s += 1;
                }
            }
        s = 1/s;
        for(i=0;i<n;++i)
            {
            one = iptr[i];
            if(one == two)
                i2ptr[i] = s;
            }
        }
    for(i=0;i<n;++i)
        for(k=0;k<col;++k)
            mptr[k*n+i] = mptr[k*n+i]*i2ptr[i];
    for(j=0;j<nl;++j)
        for(i=0;i<n;++i)
            {
            one = iptr[i];
            two = uiptr[j];
            for(k=0;k<col;++k)
                if(one == two)
                    umptr[k*nl+j] = umptr[k*nl+j] + mptr[k*n+i];
            }
    return umat;            
    }
    

SEXP getuit(SEXP x, SEXP xu, SEXP n, SEXP m, SEXP check)
    {
    int ii, jj, N, M, nProtected = 0, n_comps = 0, *iptr,ch;
    double one, two, *xptr, *xuptr;
    SEXP res, ind, names;

    N = INTEGER(n)[0];
    M = INTEGER(m)[0];
    ch = INTEGER(check)[0];

    PROTECT(ind = allocVector(INTSXP, N));
    ++nProtected;
    ++n_comps;
	
    xptr = REAL(x);
    xuptr = REAL(xu);
    iptr = INTEGER(ind);

    if(ch == 1)
        {
        for(jj = 0; jj < M; ++jj)
            {
            for(ii = 0; ii < N; ++ii) 
                {
                one = xptr[ii];
                two = xuptr[jj];
                if(one == two)
                    {
                    iptr[ii] = jj + 1;
                    }
                }
            }
        }
    if(ch == 2)
        {
        double one2,two2;
        for(jj = 0; jj < M; ++jj)
            {
            for(ii = 0; ii < N; ++ii) 
                {
                one = xptr[ii];
                one2 = xptr[ii + N];
                two = xuptr[jj];
                two2 = xuptr[jj + M];
                if((one == two) && (one2 == two2))
                    {
                    iptr[ii] = jj + 1;
                    }
                }
            }
        }

    PROTECT(res = allocVector(VECSXP, n_comps));
    ++nProtected;
    SET_VECTOR_ELT(res, 0, ind);

    PROTECT(names = allocVector(STRSXP, n_comps));
    ++nProtected;
    SET_STRING_ELT(names, 0, mkChar("ind"));
    setAttrib(res, R_NamesSymbol, names);
        
    UNPROTECT(nProtected);
    return res;
    }
	
	
SEXP mycenter2(SEXP x, SEXP n)
    {
    int ii, N;
    double mymean, mysum = 0, *xptr;
	
    xptr = REAL(x);
    N = INTEGER(n)[0];

    for(ii = 0; ii < N; ++ii)
        {
        mysum += xptr[ii];
        }

    mymean = mysum/N;
	
    for(ii = 0; ii < N; ++ii)
        {
        xptr[ii] = xptr[ii] - mymean;
        }
           
    return x;
    }
    

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
    while(up < down) 
        {            
        while(x[down] >= x[up] && up < down)      
            down--;              
        if(up != down) 
            {                    
            swapd(&x[up], &x[down]);   
            up++;    
            }
        while(x[up] <= x[down] && up < down)        
            up++;                 
        if(up != down) 
            {                    
            swapd(&x[up], &x[down]);   
            down--;   
            } 
       }       
  if(start < up)   
      quicksort_body(x, start, up-1); 
  if(end > down)  
      quicksort_body(x, down+1, end);  
  }


void quicksort(int n, double *x)
    { 
    quicksort_body(x, 0, n-1);    
    }     


/*double qcheckfun(int check, int np, double *a)*/
/*	{*/
/*	double out;*/
/*	if(check==0)*/
/*		out = (a[np-1] + a[np])/2;*/
/*	else*/
/*		out = a[np-1];*/
/*	return out;*/
/*	}    */


void retrans(double *rq, double *draws, double *redraws, int k, int m)
	{
	int i,j,jj;
	double tmp = 0.0;
	for(i=0;i<m;++i)
		{
		for(j=0;j<k;++j)
			{
			tmp = 0.0;
			for(jj=0;jj<k;++jj)
				tmp += rq[j + jj*k]*draws[jj + k*i];
			redraws[j + k*i] = tmp;
			}
		}
	}


SEXP redraw(SEXP rq, SEXP draws)
	{
	int k = ncols(rq);
	int m = ncols(draws);
	SEXP redraws;
	PROTECT(redraws = allocMatrix(REALSXP, k, m));
	retrans(REAL(rq), REAL(draws), REAL(redraws), k, m);
	UNPROTECT(1);
	return redraws;
	}


SEXP qhelp(SEXP draws, SEXP basis, SEXP iter, SEXP n, SEXP k, SEXP NP1, SEXP NP2, SEXP NP3, SEXP NP4, SEXP NP5)
    {
    int i,j,ii,jj,NC,ITER,N,np1,np2,np3,np4,np5,nProtected = 0,np11,np12,np13,np14,np15;
    SEXP out,TMP,q1,q2,q3,q4,q5,q6,names;
    
    double tmp=0.0,sum=0.0,sm=0.0;
    
    N = INTEGER(n)[0];
    NC = INTEGER(k)[0];
    ITER = INTEGER(iter)[0];
    
    PROTECT(names = allocVector(STRSXP, 6));
    ++nProtected;
        
    PROTECT(out = allocVector(VECSXP, 6));
    ++nProtected;

    PROTECT(TMP = allocVector(REALSXP, ITER));
    ++nProtected;
    
    PROTECT(q1 = allocVector(REALSXP, N));
    ++nProtected;
    
    PROTECT(q2 = allocVector(REALSXP, N));
    ++nProtected;
    
    PROTECT(q3 = allocVector(REALSXP, N));
    ++nProtected;
    
    PROTECT(q4 = allocVector(REALSXP, N));
    ++nProtected;
    
    PROTECT(q5 = allocVector(REALSXP, N));
    ++nProtected;
    
    PROTECT(q6 = allocVector(REALSXP, N));
    ++nProtected;
    
    np11 = REAL(NP1)[0];
    np12 = REAL(NP1)[0];
    np13 = REAL(NP1)[0];
    np14 = REAL(NP1)[0];
    np15 = REAL(NP1)[0];
    
    PROTECT(NP1 = coerceVector(NP1, INTSXP));
    ++nProtected;
    PROTECT(NP2 = coerceVector(NP2, INTSXP));
    ++nProtected;
    PROTECT(NP3 = coerceVector(NP3, INTSXP));
    ++nProtected;
    PROTECT(NP4 = coerceVector(NP4, INTSXP));
    ++nProtected;
    PROTECT(NP5 = coerceVector(NP5, INTSXP));
    ++nProtected;
    
    np1 = ITER - INTEGER(NP1)[0] - 1;
    np2 = ITER - INTEGER(NP2)[0] - 1;
    np3 = ITER - INTEGER(NP3)[0] - 1;
    np4 = ITER - INTEGER(NP4)[0] - 1;
    np5 = ITER - INTEGER(NP5)[0] - 1;
    
    double *bptr, *dptr, *tptr, *q1ptr,*q2ptr, *q3ptr, *q4ptr, *q5ptr, *q6ptr;
    bptr = REAL(basis);
    dptr = REAL(draws);
    tptr = REAL(TMP);
    q1ptr = REAL(q1);
    q2ptr = REAL(q2);
    q3ptr = REAL(q3);
    q4ptr = REAL(q4);
    q5ptr = REAL(q5);
    q6ptr = REAL(q6);
    
    /*Rprintf("N %d\n",N);
    Rprintf("ITER %d\n",ITER);
    Rprintf("NC %d\n",NC);*/

    const int qcheck = (np11 - floor(np11));
    
    for(i=0;i<N;++i)
        {
        sum = 0.0;
        for(ii=0;ii<ITER;++ii)
            {
            tmp = 0.0;
            for(j=0;j<NC;++j)
                {
                tmp += bptr[j*N + i]*dptr[j + NC*ii];
                }
            tptr[ii] = tmp;
            sum += tmp;
            }
            
        quicksort(ITER,tptr);      
        q1ptr[i] = sum/ITER;
              
        if((np11 - floor(np11)) == 0)    
            q2ptr[i] = (tptr[np1 - 1] + tptr[np1])/2;
        else
            q2ptr[i] = tptr[np1 - 1];
        if((np12 - floor(np12)) == 0)    
            q3ptr[i] = (tptr[np2 - 1] + tptr[np2])/2;
        else
            q3ptr[i] = tptr[np2 - 1];
        if((np13 - floor(np13)) == 0)    
            q4ptr[i] = (tptr[np3 - 1] + tptr[np3])/2;
        else
            q4ptr[i] = tptr[np3 - 1];
        if((np14 - floor(np14)) == 0)    
            q5ptr[i] = (tptr[np4 - 1] + tptr[np4])/2;
        else
            q5ptr[i] = tptr[np4 - 1];
        if((np15 - floor(np15)) == 0)    
            q6ptr[i] = (tptr[np5 - 1] + tptr[np5])/2;
        else
            q6ptr[i] = tptr[np5 - 1];

/*	q2ptr[i] = qcheckfun(qcheck, np1, tptr);*/
/*	q3ptr[i] = qcheckfun(qcheck, np2, tptr);*/
/*	q4ptr[i] = qcheckfun(qcheck, np3, tptr);*/
/*	q5ptr[i] = qcheckfun(qcheck, np4, tptr);*/
/*	q6ptr[i] = qcheckfun(qcheck, np5, tptr);*/
            
        sm += q1ptr[i];
        }
        
    double ce = sm/N;
        
    for(i=0;i<N;++i)
        {     
        q1ptr[i] = q1ptr[i] - ce;
        q2ptr[i] = q2ptr[i] - ce;
        q3ptr[i] = q3ptr[i] - ce;
        q4ptr[i] = q4ptr[i] - ce;
        q5ptr[i] = q5ptr[i] - ce;
        q6ptr[i] = q6ptr[i] - ce;
        }
        
    SET_VECTOR_ELT(out, 0, q1);
    SET_VECTOR_ELT(out, 1, q2);
    SET_VECTOR_ELT(out, 2, q3);
    SET_VECTOR_ELT(out, 3, q4);
    SET_VECTOR_ELT(out, 4, q5);
    SET_VECTOR_ELT(out, 5, q6);
    
    SET_STRING_ELT(names, 0, mkChar("mean"));
    SET_STRING_ELT(names, 1, mkChar("q1"));
    SET_STRING_ELT(names, 2, mkChar("median"));
    SET_STRING_ELT(names, 3, mkChar("q2"));
    SET_STRING_ELT(names, 4, mkChar("q11"));
    SET_STRING_ELT(names, 5, mkChar("q21"));
    
    setAttrib(out, R_NamesSymbol, names);
    
    UNPROTECT(nProtected);
    return out;
    }
    
    
SEXP cpos(SEXP p, SEXP K, SEXP pos)
    {
    int i,n,k;
    n = INTEGER(K)[0];
    k = n + 1;
    
    double tmp,*pptr,asum,xsum,ysum;
    
    pptr = REAL(p);
    
    asum = 0;
    xsum = 0;
    ysum = 0;
    
    for(i=0;i<n;++i)
        {
        tmp = pptr[i]*pptr[i + k + 1] - pptr[i + 1]*pptr[i + k];
        asum += tmp;
        xsum += (pptr[i] + pptr[i + 1])*tmp;
        ysum += (pptr[i + k] + pptr[i + k + 1])*tmp;
        }
        
    tmp = 1/(3*asum);
    REAL(pos)[0] = tmp*xsum;
    REAL(pos)[1] = tmp*ysum;
    
    return pos;
    }


//void cpos(SEXP p, SEXP K, SEXP pos)
//    {
//    int i,n,k;
//    n = INTEGER(K)[0];
//    k = n + 1;
    
//    double tmp,*pptr,asum,xsum,ysum;
    
//    pptr = REAL(p);
    
//    asum = 0;
//    xsum = 0;
//    ysum = 0;
    
//    for(i=0;i<n;++i)
//        {
//        tmp = pptr[i]*pptr[i + k + 1] - pptr[i + 1]*pptr[i + k];
//        asum += tmp;
//        xsum += (pptr[i] + pptr[i + 1])*tmp;
//        ysum += (pptr[i + k] + pptr[i + k + 1])*tmp;
//        }
        
//    tmp = 1/(3*asum);
//    REAL(pos)[0] = tmp*xsum;
//    REAL(pos)[1] = tmp*ysum;
//    }
    

SEXP  getListElement(SEXP list, char *str)
    {
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;
     
    for(i=0;i<length(list);i++)
        if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) 
            {
            elmt = VECTOR_ELT(list, i);
            break;
            }
            
    return elmt;
    }
    
    
//SEXP redmat(SEXP mat)
//    {
//    int i,n,k,nProtected = 0;
//    n = length(mat);
//    k = nrows(mat);
    
//    SEXP tmp;
//    PROTECT(tmp = allocVector(REALSXP, (n-2)));
//    ++nProtected;
    
//    SEXP dim;
//    PROTECT(dim = allocVector(INTSXP, 2));
//    ++nProtected;
    
//    for(i=1;i<k;++i)
//        REAL(tmp)[i-1] = REAL(mat)[i]; 
//    for(i=k;i<(n-1);++i)
//        REAL(tmp)[i-1] = REAL(mat)[i + 1]; 
     
//    INTEGER(dim)[0] = k - 1; 
//    INTEGER(dim)[1] = 2;   
//    setAttrib(tmp, R_DimSymbol, dim);
        
//    UNPROTECT(nProtected);
//    return tmp;
//    }


SEXP change(SEXP mat)
    {
    int i,j,n,k,T,nProtected = 0;
    n = nrows(mat);
    k = ncols(mat);
    
    T = n*k + k;
    
    double *mptr,*optr;
    
    SEXP out;
    PROTECT(out = allocVector(REALSXP, T));
    ++nProtected;
    
    mptr = REAL(mat);
    optr = REAL(out);
    
    for(i=0;i<n;++i)
        {
        for(j=0;j<k;++j)
            {
            optr[j*n + i + j] = mptr[j*n + i];            
            }
        }
    for(j=0;j<k;++j)
        {
        optr[j*n + n + j] = mptr[j*n];
        }
    
    SEXP dim;
    PROTECT(dim = allocVector(INTSXP, 2));
    ++nProtected;   

    INTEGER(dim)[0] = (n+1); 
    INTEGER(dim)[1] = k;   
    setAttrib(out, R_DimSymbol, dim);
    
    UNPROTECT(nProtected);
    return out;
    }

   
SEXP myNArem(SEXP mat)
    {
    int i,j,ok,n,k,dsum,nProtected = 0,check;
    n = nrows(mat);
    k = ncols(mat);
    
    double *mptr,*nmptr;
    int *iptr;
    
    SEXP ind;
    PROTECT(ind = allocVector(INTSXP, n));
    ++nProtected;
    
    dsum = 0,ok = 0;
    iptr = INTEGER(ind);
    mptr = REAL(mat);
    
    for(i=0;i<n;++i)
        {
        check = 0;
        for(j=0;j<k;++j)
            {
            if(ISNA(mptr[j*n + i]))
                {
                check += 1;
                }
            }
        iptr[i] = check;
        if(iptr[i] > 0)
            dsum += 1;
        }

    int T = n - dsum;
    
    SEXP newmat;
    PROTECT(newmat = allocVector(REALSXP,T*k));
    ++nProtected;
    
    nmptr = REAL(newmat);
    
    for(i=0;i<n;++i)
        {
        if(iptr[i] < 1)
            {
            for(j=0;j<k;++j)
                {
                nmptr[j*T + ok] = mptr[j*n + i];
                }
            ok += 1;
            }
        } 
    
    SEXP dim;
    PROTECT(dim = allocVector(INTSXP, 2));
    ++nProtected;   
    
    INTEGER(dim)[0] = T; 
    INTEGER(dim)[1] = k;   
    setAttrib(newmat, R_DimSymbol, dim);
        
    UNPROTECT(nProtected);
    return newmat;
    }

    
SEXP cdist(SEXP map, SEXP N, SEXP b)
    {
    int i,j,n,nProtected = 0,nn1,nn2;
    n = INTEGER(N)[0];
    
    double dx,dy;

    SEXP p1 = R_NilValue;
    PROTECT(p1);
    ++nProtected;
    
    SEXP p2 = R_NilValue;
    PROTECT(p2);
    ++nProtected;
    
    SEXP pos1;
    PROTECT(pos1 = allocVector(REALSXP, 2));
    ++nProtected;
    
    SEXP pos2;
    PROTECT(pos2 = allocVector(REALSXP, 2));
    ++nProtected;

    SEXP NP1;
    PROTECT(NP1 = allocVector(INTSXP, 1));
    ++nProtected;
    
    SEXP NP2;
    PROTECT(NP2 = allocVector(INTSXP, 1));
    ++nProtected;
    
    SEXP tmp1;
    PROTECT(tmp1 = allocVector(REALSXP, 2));
    ++nProtected;
    REAL(tmp1)[0] = 0;
    REAL(tmp1)[1] = 0;
    
    SEXP tmp2;
    PROTECT(tmp2 = allocVector(REALSXP, 2));
    ++nProtected;
    REAL(tmp2)[0] = 0;
    REAL(tmp2)[1] = 0;
    
    SEXP names;
    PROTECT(names = allocVector(STRSXP, 2));
    ++nProtected;
     
    SEXP out;   
    PROTECT(out = allocVector(VECSXP, 2));
    ++nProtected;
    
    SEXP ctr;   
    PROTECT(ctr = allocVector(REALSXP, n*2));
    ++nProtected;
    
    SEXP dimb;
    PROTECT(dimb = allocVector(INTSXP, 2));
    ++nProtected;
    
    SEXP dimctr;
    PROTECT(dimctr = allocVector(INTSXP, 2));
    ++nProtected;
    
    int *np1ptr,*np2ptr;
    double *ctrptr,*bptr;
    
    np1ptr = INTEGER(NP1);
    np2ptr = INTEGER(NP2);
    ctrptr = REAL(ctr);
    bptr = REAL(b);
    
    int ii,kk;
    double tmp,*pptr,asum,xsum,ysum;
    
    for(i=0;i<n;++i)
        {
        for(j=0;j<n;++j)
            {
            if((i!=j) && (j>i))
                {               
                SET_VECTOR_ELT(map, i, myNArem(VECTOR_ELT(map, i)));
                np1ptr[0] = nrows(VECTOR_ELT(map, i)) - 1;

                SET_VECTOR_ELT(map, j, myNArem(VECTOR_ELT(map, j)));
                np2ptr[0] = nrows(VECTOR_ELT(map, j)) - 1;
                    
                nn1 = np1ptr[0] + 1;
                nn2 = np2ptr[0] + 1;
                    
                if(REAL(VECTOR_ELT(map, i))[0] != REAL(VECTOR_ELT(map, i))[nn1 - 1] || REAL(VECTOR_ELT(map, i))[nn1] != REAL(VECTOR_ELT(map, i))[2*nn1 - 1])
                      {
                      SET_VECTOR_ELT(map, i, change(VECTOR_ELT(map, i)));
                      np1ptr[0] = nrows(VECTOR_ELT(map, i)) - 1;
                      }
                if(REAL(VECTOR_ELT(map, j))[0] != REAL(VECTOR_ELT(map, j))[nn2 - 1] || REAL(VECTOR_ELT(map, j))[nn2] != REAL(VECTOR_ELT(map, j))[2*nn2 - 1])
                      {
                      SET_VECTOR_ELT(map, j, change(VECTOR_ELT(map, j)));
                      np2ptr[0] = nrows(VECTOR_ELT(map, j)) - 1;
                      }
                
                //pos1 = cpos(p1,NP1,tmp1);
                //pos2 = cpos(p2,NP2,tmp2);
                
                asum = 0;
                xsum = 0;
                ysum = 0;
                kk = np1ptr[0] + 1;
                
                pptr = REAL(VECTOR_ELT(map, i));
    
                for(ii=0;ii<np1ptr[0];++ii)
                    {
                    tmp = pptr[ii]*pptr[ii + kk + 1] - pptr[ii + 1]*pptr[ii + kk];
                    asum += tmp;
                    xsum += (pptr[ii] + pptr[ii + 1])*tmp;
                    ysum += (pptr[ii + kk] + pptr[ii + kk + 1])*tmp;
                    }
        
                tmp = 1/(3*asum);
                REAL(pos1)[0] = tmp*xsum;
                REAL(pos1)[1] = tmp*ysum;
                
                asum = 0;
                xsum = 0;
                ysum = 0;
                kk = np2ptr[0] + 1;
               
                pptr = REAL(VECTOR_ELT(map, j));
  
                for(ii=0;ii<np2ptr[0];++ii)
                    {
                    tmp = pptr[ii]*pptr[ii + kk + 1] - pptr[ii + 1]*pptr[ii + kk];
                    asum += tmp;
                    xsum += (pptr[ii] + pptr[ii + 1])*tmp;
                    ysum += (pptr[ii + kk] + pptr[ii + kk + 1])*tmp;
                    }
        
                tmp = 1/(3*asum);
                REAL(pos2)[0] = tmp*xsum;
                REAL(pos2)[1] = tmp*ysum;
                
                ctrptr[i] = REAL(pos1)[0];
                ctrptr[n + i] = REAL(pos1)[1];
                
                if(j == (n - 1))
                    {
                    ctrptr[j] = REAL(pos2)[0];
                    ctrptr[n + j] = REAL(pos2)[1];
                    }
                
                dx = REAL(pos1)[0] - REAL(pos2)[0];
                dy = REAL(pos1)[1] - REAL(pos2)[1];
               
                bptr[j*n + i] = pythag(dx,dy);
                bptr[i*n + j] = bptr[j*n + i];
                }
            }
        }
     
    INTEGER(dimb)[0] = n; 
    INTEGER(dimb)[1] = n;   
    setAttrib(b, R_DimSymbol, dimb);
    
    INTEGER(dimctr)[0] = n;   
    INTEGER(dimctr)[1] = 2;   
    setAttrib(ctr, R_DimSymbol, dimctr);
        
    SET_VECTOR_ELT(out, 0, b);
    SET_VECTOR_ELT(out, 1, ctr);
    
    SET_STRING_ELT(names, 0, mkChar("distance"));
    SET_STRING_ELT(names, 1, mkChar("centroids"));
    
    setAttrib(out, R_NamesSymbol, names);    
        
    UNPROTECT(nProtected);
    return out;
    }


SEXP cprobs(SEXP S2, SEXP LAMBDA, SEXP BETA, SEXP MBETA, SEXP S, SEXP POS, SEXP N, SEXP K)
    {
    int i,j,d,n,pos,k,nProtected = 0;
    n = INTEGER(N)[0];
    k = INTEGER(K)[0];
    pos = n - INTEGER(POS)[0] - 1;
    double ld,p,db1,db2,tmp1,tmp2,check = 0,smp = n; 
    
    // ntmp = (n*(n+1))/2 
    
    SEXP C1;   
    PROTECT(C1 = allocVector(REALSXP, n));
    ++nProtected;
    
    SEXP C2;   
    PROTECT(C2 = allocVector(REALSXP, n));
    ++nProtected;
    
    SEXP probs;
    PROTECT(probs = allocVector(REALSXP, 1));
    ++nProtected;
    
    double *s2ptr,*laptr,*sptr,*betaptr,*mbetaptr,*c1ptr,*c2ptr;
    s2ptr = REAL(S2);
    laptr = REAL(LAMBDA);
    sptr = REAL(S);
    betaptr = REAL(BETA);
    mbetaptr = REAL(MBETA);
    c1ptr = REAL(C1);
    c2ptr = REAL(C2);
    
    for(j=0;j<n;++j)
        {
        for(i=0;i<n;++i)
            {
            if(j > i)
                {
                tmp1 = 0;
                tmp2 = 0;
                ld = 1;
            
                for(d=0;d<k;++d)
                    {
                    p = s2ptr[j]/(1.0 + laptr[j]*sptr[d]);
                    ld *= p;
                    db1 = betaptr[i*k + d] - mbetaptr[j*k + d];
                    db2 = 0.0 - mbetaptr[j*k + d];
                    tmp1 += 0.5*(db1*p*db1);
                    tmp2 += 0.5*(db2*p*db2);
                    }
                ld = log(ld);
                }
            c1ptr[i] = 0.5*ld - tmp1;
            c2ptr[i] = 0.5*ld - tmp2;
            }
      
        quicksort(n,c1ptr);
        quicksort(n,c2ptr); 
        
        if(c1ptr[pos] < c2ptr[pos])
            check += 1;
        }
    REAL(probs)[0] = check/smp;
        
    UNPROTECT(nProtected);
    return probs;
    }
    

SEXP cpoint(SEXP poly1, SEXP poly2)
    {
    int i,j,np1,np2;
    
    np1 = nrows(poly1)-1;
    np2 = nrows(poly2)-1;
    
    double *pptr1,*pptr2;
    pptr1 = REAL(poly1);
    pptr2 = REAL(poly2);
    
    SEXP out;
    PROTECT(out = allocVector(INTSXP, 1));
    INTEGER(out)[0] = 0;
    
    for(i=0;i<np1;++i)
        for(j=0;j<np2;++j)
            if(pptr1[i] == pptr2[j])
                if(pptr1[i + np1] == pptr2[j + np2])
                    {
                    INTEGER(out)[0] = 1;
                    break;
                    }
     
    UNPROTECT(1);
    return out;
    }
    

SEXP qhelp2(SEXP draws, SEXP basis, SEXP iter, SEXP n, SEXP k, 
	    SEXP NP1, SEXP NP2, SEXP NP3, SEXP NP4, SEXP NP5,
	    SEXP response, SEXP eta, SEXP ind)
    {
    int i,j,ii,jj,NC,ITER,N,nr,np1,np2,np3,np4,np5,nProtected = 0,np11,np12,np13,np14,np15;
    SEXP out,TMP,fout;
    
    double tmp=0.0,sum=0.0;
    
    N = INTEGER(n)[0];
    NC = INTEGER(k)[0];
    ITER = INTEGER(iter)[0];
    nr = length(ind);
        
    PROTECT(out = allocMatrix(REALSXP, N, 6));
    ++nProtected;

    PROTECT(fout = allocMatrix(REALSXP, N, 9));
    ++nProtected;

    PROTECT(TMP = allocVector(REALSXP, ITER));
    ++nProtected;
    
    np11 = REAL(NP1)[0];
    np12 = REAL(NP1)[0];
    np13 = REAL(NP1)[0];
    np14 = REAL(NP1)[0];
    np15 = REAL(NP1)[0];
    
    PROTECT(NP1 = coerceVector(NP1, INTSXP));
    ++nProtected;
    PROTECT(NP2 = coerceVector(NP2, INTSXP));
    ++nProtected;
    PROTECT(NP3 = coerceVector(NP3, INTSXP));
    ++nProtected;
    PROTECT(NP4 = coerceVector(NP4, INTSXP));
    ++nProtected;
    PROTECT(NP5 = coerceVector(NP5, INTSXP));
    ++nProtected;
    
    np1 = ITER - INTEGER(NP1)[0] - 1;
    np2 = ITER - INTEGER(NP2)[0] - 1;
    np3 = ITER - INTEGER(NP3)[0] - 1;
    np4 = ITER - INTEGER(NP4)[0] - 1;
    np5 = ITER - INTEGER(NP5)[0] - 1;
    
    double *bptr,*dptr,*tptr,*outptr,*responseptr,*etaptr,*foutptr;
    int *indptr;
    bptr = REAL(basis);
    dptr = REAL(draws);
    tptr = REAL(TMP);
    outptr = REAL(out);
    responseptr = REAL(response);
    etaptr = REAL(eta);
    foutptr = REAL(fout);
    indptr = INTEGER(ind);

    const int qcheck = (np11 - floor(np11));
    
    for(i=0;i<N;++i)
        {
        sum = 0.0;
        for(ii=0;ii<ITER;++ii)
            {
            tmp = 0.0;
            for(j=0;j<NC;++j)
                {
                tmp += bptr[j*N + i]*dptr[j + NC*ii];
                }
            tptr[ii] = tmp;
            sum += tmp;
            }
            
        quicksort(ITER,tptr);      
        outptr[i] = sum/ITER;
              
        if((np11 - floor(np11)) == 0)    
            outptr[i + 5*N] = (tptr[np1 - 1] + tptr[np1])/2;
        else
            outptr[i + 5*N] = tptr[np1 - 1];
        if((np12 - floor(np12)) == 0)    
            outptr[i + 3*N] = (tptr[np2 - 1] + tptr[np2])/2;
        else
            outptr[i + 3*N] = tptr[np2 - 1];
        if((np13 - floor(np13)) == 0)    
            outptr[i + 1*N] = (tptr[np3 - 1] + tptr[np3])/2;
        else
            outptr[i + 1*N] = tptr[np3 - 1];
        if((np14 - floor(np14)) == 0)    
            outptr[i + 4*N] = (tptr[np4 - 1] + tptr[np4])/2;
        else
            outptr[i + 4*N] = tptr[np4 - 1];
        if((np15 - floor(np15)) == 0)    
            outptr[i + 2*N] = (tptr[np5 - 1] + tptr[np5])/2;
        else
            outptr[i + 2*N] = tptr[np5 - 1];
        }
        
    int get = 0;

    for(i=0;i<nr;++i)
        {
	get = indptr[i] - 1;
	for(j=0;j<6;++j)
		{
		foutptr[i + j*nr] = outptr[get + j*N];
		}
	foutptr[i + 7*N] = (foutptr[i + 1*N] < 0 && foutptr[i + 5*N] < 0)*(-1) + (foutptr[i + 1*N] <= 0 & foutptr[i + 5*N] >= 0)*0 + (foutptr[i + 1*N] > 0 && foutptr[i + 5*N] > 0)*1;
	foutptr[i + 8*N] = (foutptr[i + 2*N] < 0 && foutptr[i + 4*N] < 0)*(-1) + (foutptr[i + 2*N] <= 0 && foutptr[i + 4*N] >= 0)*0 + (foutptr[i + 2*N] > 0 && foutptr[i + 4*N] > 0)*1; 
	foutptr[i + 6*N] = responseptr[i] - etaptr[i] + foutptr[i]; 
        }
    
    UNPROTECT(nProtected);
    return fout;
    }
