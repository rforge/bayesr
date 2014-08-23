/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */


// Modul zur Erzeugung von Zufallszahlen



#include "Random.h"


// BEGIN: DSB //

// define what is NAN, and how we determine whether a double is infinite.
#if defined(MICROSOFT_VISUAL)
#include <limits>
    bool
    infinite(double x)
    {
        return ABS(x) > DBL_MAX;
    }
#elif defined(__BUILDING_GNU)
    bool
    infinite(double x)
    {
        return ABS(x) > DBL_MAX;
    }
#else
#include"../values.h"
    bool
    infinite(double x)
    {
        return ABS(x) > DBL_MAX;
    }
#endif

#ifndef NAN
  #define NAN DBL_MAX
#endif

// define a "repeat" syntax (for the new generators adapted from R)
#define repeat for(;;)

// END: DSB //




namespace randnumbers
{

double Phi(const double & x)
{
if (x==0)
    return 0.5;
  else
    {
    double a=0.0;
    double b=0.0;
    if (x > 0)
      {
      a = 0;
      b = x;
      }
    else if (x<0)
      {
      a=x;
      b=0;
      }

    double h = (b-a)/50.0;
    double x2 = -0.5*a*a;
    double sum= exp(x2);
    double xhelp = a+h;
    unsigned i;
    for (i=1;i<=25;i++)
      {
      x2 = -0.5*xhelp*xhelp;
      sum+=4*exp(x2);
      xhelp+=h;

      x2 = -0.5*xhelp*xhelp;
      if (i == 25)
        sum+=exp(x2);
      else
        sum+=2*exp(x2);
      xhelp+=h;

      }

    double c = 0.13298076*h;

    sum *= c;

    if (x > 0)
      return sum+0.5;
    else
      return 0.5-sum;
    }

}


double Phi2(const double & x)
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


double __EXPORT_TYPE invPhi (const double & p);


double invPhi2 (const double & p)
  {

    double pt;
    if (p < 0.5)
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

    if (p<0.5)
      return -t+r1/r2;
    else
      return t-r1/r2;

  }


double uniform(void)
  {


/*
  randomize();
  static double  m = 2147483648;
  static double a = 843314861;
  static double c = 453816693;


  static double x = int(rand()+5);

  x = fmod(a*x+c,m);
  return x/m;
*/

  int zufall = 0;
  while ((zufall == 0) || (zufall == RAND_MAX))
	 {
	 zufall = rand();
	 }
  return double(zufall)/double(RAND_MAX);

  }

double uniform_ab(double a, double b)
{
    double zufall = 0;
    zufall = a+uniform()*(b-a);
    return(zufall);
}



double rand_gamma(double a,double b)
  {
  if (a > 1)
	 {
	 double h1 = a-1;         // h1 entspricht b in Devroye (1986)
	 double h2 = 3*a-0.75;    // h2 entspricht c in Devroye (1986)
	 double U,V,W,Y,X,Z;
	 int accept = 0;
	 do
		{
		U = uniform();
		V = uniform();
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
		double X = rand_gamma(a+1,1)*pow(uniform(),1/a);
		return X/b;
		}
	 }
  }


/*
double Phi(const double & x)
  {
  if (x==0)
    return 0.5;
  else
    {
    double a,b;
    if (x > 0)
      {
      a = 0;
      b = x;
      }
    else if (x<0)
      {
      a=x;
      b=0;
      }

    double h = (b-a)/50.0;
    double x2 = -0.5*a*a;
    double sum= exp(x2);
    double xhelp = a+h;
    unsigned i;
    for (i=1;i<=25;i++)
      {
      x2 = -0.5*xhelp*xhelp;
      sum+=4*exp(x2);
      xhelp+=h;

      x2 = -0.5*xhelp*xhelp;
      if (i == 25)
        sum+=exp(x2);
      else
        sum+=2*exp(x2);
      xhelp+=h;

      }

    double c = 0.13298076*h;

    sum *= c;

    if (x > 0)
      return sum+0.5;
    else
      return 0.5-sum;
    }

  }

*/

double ksdist(int kmax, double lambda)
  {
  double result = 0.0;
  int k;
  for(k=-kmax; k<=kmax; k++)
    if((k%2)==0)
      result += exp(-2*SQR(lambda)*SQR(k));
    else
      result -= exp(-2*SQR(lambda)*SQR(k));
  return(result);
}

double kssample(void)
  {
  double U = uniform();
  int accept = 0;
  double E0, E1, n, G, X, W, Z, P, Q, U2, E, U3, dummy;
  int j;
  if(U < Ft)  // generator for the leftmost interval
    {
    do{

      do{
        E0 = -log(uniform());        // exponential
        E1 = -log(uniform());       // exponential
        E0 = E0/(1-1/(2*ts));
        E1 = 2*E1;
        G = ts+E0;
        accept = (SQR(E0) <= (ts*E1*(G+ts)));

        if(!accept)
          accept = ((G/ts-1-log(G/ts)) <= (E1));
        }
      while (!accept);

      X = PI/sqrt(8*G);
      W = 0.0;
      Z = 1/(2*G);
      P = exp(-G);
      n = 1.0;
      Q = 1.0;
      U2 = uniform();

      do{
        W = W + Z*Q;
        if(U2 >= W)
          return(X);
        n = n+2;
        Q = P;
      for(j=2; j<=(SQR(n)-1); j++)
         Q *= P;
         W = W - SQR(n)*Q;
         }
      while(U2 >= W);

      }
    while(1);
    }
  else  // generator for the rightmost interval
    {

    do
      {
      E = -log(uniform()); // exponential
      U3 = uniform();
      X = sqrt(SQR(t) + E/2.0);
      W = 0.0;
      n = 1.0;
      Z = exp(-2*SQR(X));
      do
        {
        n++;
        dummy = Z;
        for(j=2; j<=(SQR(n)-1); j++)
          dummy *= Z;
        W = W + SQR(n)*dummy;
        if (U3 >= W)
          return(X);
          n++;
        dummy = Z;
        for(j=2; j<=(SQR(n)-1); j++)
          dummy *= Z;
        W = W - SQR(n)*dummy;

        }
      while(U3 > W);
      }
    while(1);
    }

  }

double rand_normal(void)
  {
  double u1 = uniform();
  double u2 = uniform();
  return sqrt(-2*log(u1))*sin(6.2831853*u2);
  }

double trunc_normal(const double & a,const double & b,const double & mu,
                   const double & s)
  {
  bool accept = false;
  double rand;
  while (!accept)
    {
    rand = mu+s*rand_normal();
    if ((rand <= b) && (rand >= a))
      accept = true;
    }

  return rand;
  }

double trunc_normal2(const double & a,const double & b,const double & mu,
                    const double & s)
  {
  double at = Phi2((a-mu)/s);
  double bt = Phi2((b-mu)/s);
  double u = at+(bt-at)*uniform();
  double r = mu+s*invPhi2(u);
  if (r < a)
    r = a+0.00000001;
  if (r > b)
    r = b-0.00000001;

  return r;
  }

double trunc_normal3(const double & a,const double & b,const double & mu,
                    const double & s)
  {
  double z;
  double at = (a-mu)/s;
  double bt = (b-mu)/s;
  double u = 1;
  double r = 0;
  while(u > r)
    {
    z = (bt-at)*uniform()+at;

    if(0 < at)
      r = exp((at*at-z*z)/2);
    else if (0 > bt)
      r = exp((bt*bt-z*z)/2);
    else
      r = exp((-z*z)/2);

    u = uniform();
    }

  return mu+s*z;
  }

double trunc_normal4(const double & a, const double & mu, const double & s)
  {
  double at = Phi2((a-mu)/s);
  double u = at+(1-at)*uniform();
  double r = mu+s*invPhi2(u);
  if (r < a)
    r = a+0.00000001;

  return r;
  }

double truncnormal(const double & a,const double & b)
  {
  bool accept = false;
  double rand;
  if (a > 2.5)
    {
    double u;
    while (!accept)
      {
      u = uniform();
      rand = a+(b-a)*uniform();
      if ( u <= (Phi(rand)/Phi(a)) )
        accept=true;
      }
    }
  else if (b < -2.5)
    {
    double u;
    while (!accept)
      {
      u = uniform();
      rand = a+(b-a)*uniform();
      if ( u <= (Phi(rand)/Phi(b)) )
        accept=true;
      }
    }
  else
    {

    while (!accept)
      {
      rand = rand_normal();
      if ((rand <= b) && (rand >= a))
        accept = true;
      }

    }

  return rand;
  }



// Erzeugen eines standardnormalverteilten Zufallsvektors (Spaltenvektor)
// mit Dimension dim !

Matrix<double> rand_normvek(unsigned dim)
  {
  Matrix<double> stnorm(dim,1);
  for (unsigned i=0; i < dim; i++)
	 stnorm(i,0) = rand_normal();
  return stnorm;
  }

// Erzeugung einer Wishart verteilten Zufallsmatrix mit n Freiheitsgraden
// und Skalenparameter Sigma und Dimension q x q
// w gibt an, ob es sich bei Sigma schon um die Choleskyzerlegung einer
// Kovarianzmatrix handelt, oder ob Sigma noch zerlegt werden muss
// w = 1 (default) entspricht schon zerlegt ansonsten unzerlegt

void rand_wishart(Matrix<double> & Sigma,const unsigned & n,Matrix<double> & res)
  {
  unsigned p = Sigma.rows();
  Matrix<double> V(p,p);
  Matrix<double> E(p,p);

  unsigned i,j,k;

  double zeta;
  double sum;

  for(i=0;i<p;i++)
    for(j=0;j<p;j++)
      E(i,j) = rand_normal();

  for(i=0;i<p;i++)
    {

    zeta = rand_chisquare(n-i);

    sum=0;
    for(k=0;k<i;k++)
      sum += E(k,i)*E(k,i);

    V(i,i) = zeta + sum;

    for(j=i+1;j<p;j++)
      {
      sum =  0;
      for(k=0;k<i;k++)
        sum += E(k,i)*E(k,j);
      V(i,j) = E(i,j)*sqrt(zeta) + sum;
      V(j,i) = V(i,j);
      }

    }

  Matrix<double> R = Sigma.root();

  res = R*V*R.transposed();

  }


  //Inverse Gaussian random numbers: Devroye 1986
double rand_inv_gaussian(const double mu, const double lambda)
{
    double N;
    double Y;
    double X1;
    N = rand_normal();
    Y = N*N;
    X1 = mu + mu*mu*Y/(2*lambda)-mu*sqrt(4*mu*lambda*Y+mu*mu*Y*Y)/(2*lambda);
    if(uniform()<= mu/(mu+X1))
    {
        return X1;
    }
    else
    {
        return mu*mu/X1;
    }
}


double rand_variance(const double f)
  {
  double length = f - 1/f;
  if (f == 1.0)
    return 1.0;
  if (uniform() < length/(length+2*log(f)))
    return (1/f + length*uniform());
  else
    return pow(f, 2.0*uniform()-1.0);
  }


double invlogit(double x)
  {
  return 1/(1+exp(-x));
  }


double logit(double x)
  {
  return log(x/(1-x));
  }


double trunc_logistic_left(double mean)
  {

  double Fa, arg, result;

  Fa = invlogit(mean);
  arg = Fa + uniform()*(1.0-Fa);
  result = mean - logit(arg);
  return result;
  }


double trunc_logistic(double mean, int left)
  {
  assert(left == 1 || left == 0);

  if(left == 1)
    return trunc_logistic_left(mean);
  else
    return -trunc_logistic_left(-mean);
  }


double IG(double mu, double lambda)
  {

  double N, Y, X1, mu2;

  mu2 = SQR(mu);
  N = rand_normal();
  Y = SQR(N);
  X1 = mu + mu2*Y/(2*lambda)-mu/(2*lambda)*
       sqrt(4*mu*lambda*Y+mu2*SQR(Y));
  if(uniform() <= mu/(mu+X1))
    return X1;
  else
    return mu2/X1;
  }


double GIG(double chi)
  {
  // returns a sample from a generalized inverse gaussian
  // distribution with parameters lambda=0.5, chi, psi=1
  // Devroye, page 478 parametrization
  // assume lambda > 0

  double r = sqrt(chi);
  return r/IG(1., r);

  }


double GIG(double lambda, double psi, double chi)
  {

  double out;

  if (chi == 0)
    return rand_gamma(lambda,psi/2);
  else if (psi == 0)
    return rand_invgamma(-lambda,chi/2);
  else
    {

    double h = lambda;
    double b = sqrt(chi * psi);
    double m = ( h-1+sqrt( pow(h-1,2)+ pow(b,2)) ) / b;
    double m2 = pow(m,2);
    double log_1_over_pm = -(h-1)/2*log(m) + b/4*(m + (1/m));
    double r = (6*m + 2*h*m - b* m2 + b)/(4* m2);
    double s = (1 + h - b*m)/(2*m2);
    double p = (3*s - pow(r,2))/3;
    double q = (2* pow(r,3))/27 - (r*s)/27 + b/(-4*m2);
    double eta = sqrt(-(pow(p,3))/27);
    double y1  = 2*exp(log(eta)/3) * cos(acos(-q/(2*eta))/3) - r/3;
    double y2  = 2*exp(log(eta)/3) * cos(acos(-q/(2*eta))/3 + 2/3* 3.1415) - r/3;

    double ym;

    if ( ((h<=1.0) && (b<=1.0)) || (fabs(q/eta)>2.0) || (y1<0.0) || (y2>0.0) )
      {
      ym = (-h-1 + sqrt( pow(h+1,2) + pow(b,2)))/b;

      double a = exp(-0.5*h*log(m*ym) + 0.5*log(m/ym) + b/4*(m + 1/m - ym - 1.0/ym));

      double u = uniform();
      double v = uniform();
      double out = a * (v/u);
      while ( log(u) > (h-1)/2*log(out) - b/4*(out + 1/out) + log_1_over_pm  )
        {
        u = uniform();
        v = uniform();
        out = a * (v/u);
        }
      }
    else
      {

      double vplus = exp(log_1_over_pm + log(1/y1) + (h-1)/2*log(1/y1 + m) -
                           b/4*(1/y1 + m + 1/(1/y1 + m)) );
      double vminus = -exp(log_1_over_pm + log(-1/y2) + (h-1)/2*log(1/y2 + m) -
                           b/4*(1/y2 + m + 1/(1/y2 + m)) );

      double u = uniform();
      double v = vminus + (vplus - vminus) * uniform();
      double z = v/u;

      while (z < -m )
        {
        u = uniform();
        v = vminus + (vplus - vminus) * uniform();
        z = v/u;
        }
      out = z + m;

      while (log(u) > (log_1_over_pm + (h-1)/2*log(out) - b/4*(out + 1/out)) )
        {
        u = uniform();
        v = vminus + (vplus - vminus) * uniform();
        z = v /u;

        while (z < -m)
          {
          u = uniform();
          v = vminus + (vplus - vminus) * uniform();
          z = v/u;
          }

          out = z + m;

        }

       }

     }

  return sqrt( chi / psi ) * out;

  }


double f1old(double x, int j)
  {
  // first series approximation, guranteed to be monotone for x > 1.333
  double a;

  a = SQR(j+1)*exp(-0.5*x*(SQR(j+1)-1));
  return a;
  }


double f2old(double x, int j)
  {

  // second series approximation, guranteed to be monotone for x < 1.333
  double a;

  if ((j%2)==1)
    // odd
   a = (x/PI2)*exp(-((SQR(j)-1.)*PI2)/(2*x));
  else
    // even
   a =  SQR(j+1)*exp(-((SQR(j+1)-1.)*PI2)/(2*x));

  return a;
  }


double lambda_fc(double chi)
  {

  double fn_val,u;
  int ok_out;

  ok_out = 0;


  int  j, t=0;
  double aa, upper, lower, factor;

   while (!ok_out){

     fn_val = GIG(chi);

     // u is the random variable used to test for acceptance prob, 4 is the
     // supremum under the generalised-inverse-gaussian density
     u = uniform()*4.0;

     // using Devroye (and transformation of RVs) we can find an alternating
     // series approximation. In particular we can find a series of monotone
     // "squeezing" functions that bound the true ratio
     // (true_function(x)/sampling_function(x)) from above and below ever
     // more closely. However as in Devroye, in order to do this we need to
     // split up the x space into two regions within which we find a guaranteed
     // monotone series.

     // so in the first region.....

     if (fn_val > 1.334){

       j=1;
       factor = 4.0;
       upper = 1.0;
       // apply squeezing
       while (1){
         // first adjust the lower bound using the alternating series
         // f1() - given below
         lower = upper - f1old(fn_val,j);
         if (u < factor*lower){
           // if the draw of the acceptance prob is below the lower bound
           // then you are definatly a draw from the density - ACCEPT
           ok_out = 1;
           break;
         }
         // now adjust the upper bound
         upper = lower + f1old(fn_val,j+1);
         if (u > factor*upper){
           // if the draw of the acceptance prob is above the upper bound
           // then you are definatly NOT a draw from the density - REJECT
           ok_out=0;
           break;
         }
         // else u lies somewhere inbetween (lower, upper) so we're not sure
         // and we must continue
         j+=2;
       }
     }
     else{ // you are in the other region (the above series f1() is not
           // guarenteed to be monotone hence we must find other monotone series
       j=1;

       // aa is simply the supremum under the rejection sampler which lies at
       // the boundary
       aa = 0.5*log(2*PI)+2*log(PI)+log(4.0)-2.5*log(fn_val)-
           (SQR(PI)/(2*fn_val))+0.5*fn_val;
       factor = exp(aa);
       upper = 1;
       while (1){
         // this bit is the same as above but we use the series f2()
         lower = upper - f2old(fn_val,j);
         if (u < factor*lower){
           ok_out = 1;
           break;
         }
         upper = lower + f2old(fn_val,j+1);
         if (u > factor*upper){
           ok_out=0;
           break;
         }
         j=j+2;
       }
     }

     // t counts the number of reject steps in the sampler
     t += 1;
   }
   return fn_val;
}

double rand_beta(double a, double b)
  {
  double randgamma_a;
  double randgamma_b;
  double randbeta;

  randgamma_a = rand_gamma(a,1);
  randgamma_b = rand_gamma(b,1);

  randbeta = randgamma_a/(randgamma_a+randgamma_b);
  return randbeta;
  }

vector<double> rand_dirichlet(double nrpar, vector<double> alpha)
  {
  vector<double> randgamma(nrpar);
  vector<double> randdirichlet(nrpar);
  double randgammasum = 0;
  double randdirichletsum = 0;

  for(int i=0; i<nrpar-1; i++)
    {
     randgamma[i] = rand_gamma(alpha[i],1);
     randgammasum += randgamma[i];
    }
  for(int i=0; i<nrpar-2; i++)
    {
    randdirichlet[i] = randgamma[i]/randgammasum;
    randdirichletsum = randdirichletsum + randdirichlet[i];
    }
   randdirichlet[nrpar-1] = 1-randdirichletsum;

   return randdirichlet;
  }

unsigned bernoulli(double & prob)
  {
  unsigned res = 0;
  double u = uniform();
  if(u<=prob)
    res=1;

  return res;
  }

// BEGIN: DSB //

// adapted from R-2.10.1/src/nmath/rbinom.c

double rand_binom(double nin, double prob)
{
    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;

    static double psave = -1.0;
    static int nsave = -1;
    static int m;

    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i,ix,k, n;

    if(infinite(nin))
        return NAN;

    r = floor(nin + 0.5);
    if (r != nin)
        return NAN;

    if(infinite(prob) ||
        /* n=0, p=0, p=1 are not errors <TSL>*/
        r < 0 || prob < 0. || prob > 1.)
        return NAN;

    if (r == 0 || prob == 0.)
        return 0;
    if (prob == 1.)
        return r;

    if (r >= INT_MAX)/* evade integer overflow,
                                and r == INT_MAX gave only even values */
        return NAN; // qbinom(unif_rand(), r, pp, /*lower_tail*/ 0, /*log_p*/ 0);
    /* else */
    n = r;

//    p = fmin(prob, 1. - prob);
    p = MIN(prob, 1. - prob);
//    p = prob;
//    if((1. - prob) < prob)
//      p = 1. - prob;

    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);

    /* Setup, perform only when parameters change [using static (globals): */

    /* FIXING: Want this thread safe
       -- use as little (thread globals) as possible
    */
    if (prob != psave || n != nsave) {
        psave = prob;
        nsave = n;
        if (np < 30.0) {
            /* inverse cdf logic for mean less than 30 */
            qn = pow(q, (double) n);
            goto L_np_small;
        } else {
            ffm = np + p;
            m = ffm;
            fm = m;
            npq = np * q;
            p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
            xm = fm + 0.5;
            xl = xm - p1;
            xr = xm + p1;
            c = 0.134 + 20.5 / (15.3 + fm);
            al = (ffm - xl) / (ffm - xl * p);
            xll = al * (1.0 + 0.5 * al);
            al = (xr - ffm) / (xr * q);
            xlr = al * (1.0 + 0.5 * al);
            p2 = p1 * (1.0 + c + c);
            p3 = p2 + c / xll;
            p4 = p3 + c / xlr;
        }
    } else if (n == nsave) {
        if (np < 30.0)
            goto L_np_small;
    }

    /*-------------------------- np = n*p >= 30 : ------------------- */

    repeat {
      u = uniform() * p4;
      v = uniform();
      /* triangular region */
      if (u <= p1) {
          ix = xm - p1 * v + u;
          goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
          x = xl + (u - p1) / c;
          v = v * c + 1.0 - ABS(xm - x) / p1;
          if (v > 1.0 || v <= 0.)
              continue;
          ix = x;
      } else {
          if (u > p3) { /* right tail */
              ix = xr - log(v) / xlr;
              if (ix > n)
                  continue;
              v = v * (u - p3) * xlr;
          } else {/* left tail */
              ix = xl + log(v) / xll;
              if (ix < 0)
                  continue;
              v = v * (u - p2) * xll;
          }
      }
      /* determine appropriate way to perform accept/reject test */
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
          /* explicit evaluation */
          f = 1.0;
          if (m < ix) {
              for (i = m + 1; i <= ix; i++)
                  f *= (g / i - r);
          } else if (m != ix) {
              for (i = ix + 1; i <= m; i++)
                  f /= (g / i - r);
          }
          if (v <= f)
              goto finis;
      } else {
          /* squeezing using upper and lower bounds on log(f(x)) */
          amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
          ynorm = -k * k / (2.0 * npq);
          alv = log(v);
          if (alv < ynorm - amaxp)
              goto finis;
          if (alv <= ynorm + amaxp) {
              /* stirling's formula to machine accuracy */
              /* for the final acceptance/rejection test */
              x1 = ix + 1;
              f1 = fm + 1.0;
              z = n + 1 - fm;
              w = n - ix + 1.0;
              z2 = z * z;
              x2 = x1 * x1;
              f2 = f1 * f1;
              w2 = w * w;
              if (alv <=
                      xm * log(f1 / x1) +
                      (n - m + 0.5) * log(z / w) +
                      (ix - m) * log(w * p / (x1 * q)) +
                      (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 +
                      (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 +
                      (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 +
                      (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.0
                  )
                  goto finis;
          }
      }
  }

 L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */

  repeat {
     ix = 0;
     f = qn;
     u = uniform();
     repeat {
         if (u < f)
             goto finis;
         if (ix > 110)
             break;
         u -= f;
         ix++;
         f *= (g / ix - r);
     }
  }

 finis:
    if (psave > 0.5)
         ix = n - ix;
  return (double)ix;
}


// adapted from R-2.10.1/src/nmath/rpois.c

double rand_pois(double mu)
{
    // use static constants instead of preprocessor macros.
    static const double a0 = -0.5;
    static const double a1 = 0.3333333;
    static const double a2 = -0.2500068;
    static const double a3 = 0.2000118;
    static const double a4 = -0.1661269;
    static const double a5 = 0.1421878;
    static const double a6 = -0.1384794;
    static const double a7 = 0.1250060;

    static const double one_7 = 0.1428571428571428571;
    static const double one_12 = 0.0833333333333333333;
    static const double one_24 = 0.0416666666666666667;

    static const double M_1_SQRT_2PI = 0.398942280401432677939946059934;

    /* Factorial Table (0:9)! */
    const static double fact[10] =
    {
        1., 1., 2., 6., 24., 120., 720., 5040., 40320., 362880.
    };

    /* These are static --- persistent between calls for same mu : */
    static int l, m;

    static double b1, b2, c, c0, c1, c2, c3;
    static double pp[36], p0, p, q, s, d, omega;
    static double big_l;/* integer "w/o overflow" */
    static double muprev = 0., muprev2 = 0.;/*, muold    = 0.*/

    /* Local Vars  [initialize some for -Wall]: */
    double del, difmuk= 0., E= 0., fk= 0., fx, fy, g, px, py, t, u= 0., v, x;
    double pois = -1.;
    int k, kflag, big_mu, new_big_mu = 0;

    if(infinite(mu) || mu < 0)
//    if (finite(mu) != 0 || mu < 0)
        return NAN;

    if (mu <= 0.)
        return 0.;

    big_mu = mu >= 10.;
    if(big_mu)
        new_big_mu = 0;

    if (!(big_mu && mu == muprev)) {/* maybe compute new persistent par.s */

        if (big_mu) {
            new_big_mu = 1;
            /* Case A. (recalculation of s,d,l  because mu has changed):
             * The poisson probabilities pk exceed the discrete normal
             * probabilities fk whenever k >= m(mu).
             */
            muprev = mu;
            s = sqrt(mu);
            d = 6. * mu * mu;
            big_l = floor(mu - 1.1484);
            /* = an upper bound to m(mu) for all mu >= 10.*/
        }
        else { /* Small mu ( < 10) -- not using normal approx. */

            /* Case B. (start new table and calculate p0 if necessary) */

            /*muprev = 0.;-* such that next time, mu != muprev ..*/
            if (mu != muprev) {
                muprev = mu;
                m = MAX(1, (int) mu);
                l = 0; /* pp[] is already ok up to pp[l] */
                q = p0 = p = exp(-mu);
            }

            repeat {
                /* Step U. uniform sample for inversion method */
                u = uniform();
                if (u <= p0)
                    return 0.;

                /* Step T. table comparison until the end pp[l] of the
                   pp-table of cumulative poisson probabilities
                   (0.458 > ~= pp[9](= 0.45792971447) for mu=10 ) */
                if (l != 0) {
                    for (k = (u <= 0.458) ? 1 : MIN(l, m);  k <= l; k++)
                        if (u <= pp[k])
                            return (double)k;
                    if (l == 35) /* u > pp[35] */
                        continue;
                }
                /* Step C. creation of new poisson
                   probabilities p[l..] and their cumulatives q =: pp[k] */
                l++;
                for (k = l; k <= 35; k++) {
                    p *= mu / k;
                    q += p;
                    pp[k] = q;
                    if (u <= q) {
                        l = k;
                        return (double)k;
                    }
                }
                l = 35;
            } /* end(repeat) */
        }/* mu < 10 */

    } /* end {initialize persistent vars} */

/* Only if mu >= 10 : ----------------------- */

    /* Step N. normal sample */
    g = mu + s * rand_normal();/* rand_normal() ~ N(0,1), standard normal */

    if (g >= 0.) {
        pois = floor(g);
        /* Step I. immediate acceptance if pois is large enough */
        if (pois >= big_l)
            return pois;
        /* Step S. squeeze acceptance */
        fk = pois;
        difmuk = mu - fk;
        u = uniform(); /* ~ U(0,1) - sample */
        if (d * u >= difmuk * difmuk * difmuk)
            return pois;
    }

    /* Step P. preparations for steps Q and H.
       (recalculations of parameters if necessary) */

    if (new_big_mu || mu != muprev2) {
        /* Careful! muprev2 is not always == muprev
           because one might have exited in step I or S
           */
        muprev2 = mu;
        omega = M_1_SQRT_2PI / s;
        /* The quantities b1, b2, c3, c2, c1, c0 are for the Hermite
         * approximations to the discrete normal probabilities fk. */

        b1 = one_24 / mu;
        b2 = 0.3 * b1 * b1;
        c3 = one_7 * b1 * b2;
        c2 = b2 - 15. * c3;
        c1 = b1 - 6. * b2 + 45. * c3;
        c0 = 1. - b1 + 3. * b2 - 15. * c3;
        c = 0.1069 / mu; /* guarantees majorization by the 'hat'-function. */
    }

    if (g >= 0.) {
        /* 'Subroutine' F is called (kflag=0 for correct return) */
        kflag = 0;
        goto Step_F;
    }


    repeat {
        /* Step E. Exponential Sample */

        E = rand_expo(1); /* ~ Exp(1) (standard exponential) */

        /*  sample t from the laplace 'hat'
            (if t <= -0.6744 then pk < fk for all mu >= 10.) */
        u = 2 * uniform() - 1.;
        t = 1.8 + FSIGN(E, u);
        if (t > -0.6744) {
            pois = floor(mu + s * t);
            fk = pois;
            difmuk = mu - fk;

            /* 'subroutine' F is called (kflag=1 for correct return) */
            kflag = 1;

          Step_F: /* 'subroutine' F : calculation of px,py,fx,fy. */

            if (pois < 10) { /* use factorials from table fact[] */
                px = -mu;
                py = pow(mu, pois) / fact[(int)pois];
            }
            else {
                /* Case pois >= 10 uses polynomial approximation
                   a0-a7 for accuracy when advisable */
                del = one_12 / fk;
                del = del * (1. - 4.8 * del * del);
                v = difmuk / fk;
                if (ABS(v) <= 0.25)
                    px = fk * v * v * (((((((a7 * v + a6) * v + a5) * v + a4) *
                                          v + a3) * v + a2) * v + a1) * v + a0)
                        - del;
                else /* |v| > 1/4 */
                    px = fk * log(1. + v) - difmuk - del;
                py = M_1_SQRT_2PI / sqrt(fk);
            }
            x = (0.5 - difmuk) / s;
            x *= x;/* x^2 */
            fx = -0.5 * x;
            fy = omega * (((c3 * x + c2) * x + c1) * x + c0);
            if (kflag > 0) {
                /* Step H. Hat acceptance (E is repeated on rejection) */
                if (c * ABS(u) <= py * exp(px + E) - fy * exp(fx + E))
                    break;
            } else
                /* Step Q. Quotient acceptance (rare case) */
                if (fy - u * fy <= py * exp(px - fx))
                    break;
        }/* t > -.67.. */
    }
    return pois;
}

// END: DSB //


double digamma_exact ( double & x)
{
  double c = 8.5;
  double d1 = -0.5772156649;
  double r;
  double s = 0.00001;
  double s3 = 0.08333333333;
  double s4 = 0.0083333333333;
  double s5 = 0.003968253968;
  double value;
  double y;

//
//  Initialize.
//
  y = x;
  value = 0.0;
//
//  Use approximation if argument <= S.
//
  if ( y <= s )
  {
    value = d1 - 1.0 / y;
    return value;
  }
//
//  Reduce to DIGAMA(X + N) where (X + N) >= C.
//
  while ( y < c )
  {
    value = value - 1.0 / y;
    y = y + 1.0;
  }
//
//  Use Stirling's (actually de Moivre's) expansion if argument > C.
//
  r = 1.0 / y;
  value = value + log ( y ) - 0.5 * r;
  r = r * r;
  value = value - r * ( s3 - r * ( s4 - r * s5 ) );

  return value;
}


double trigamma_exact (double & x)
  {
  double a = 0.0001;
  double b = 5.0;
  double b2 =  0.1666666667;
  double b4 = -0.03333333333;
  double b6 =  0.02380952381;
  double b8 = -0.03333333333;
  double value;
  double y;
  double z;

  z = x;
//
//  Use small value approximation if X <= A.
//
  if ( x <= a )
  {
    value = 1.0 / x / x;
    return value;
  }
//
//  Increase argument to ( X + I ) >= B.
//
  value = 0.0;

  while ( z < b )
  {
    value = value + 1.0 / z / z;
    z = z + 1.0;
  }
//
//  Apply asymptotic formula if argument is B or greater.
//
  y = 1.0 / z / z;

  value = value + 0.5 *
      y + ( 1.0
    + y * ( b2
    + y * ( b4
    + y * ( b6
    + y *   b8 )))) / z;

  return value;
}


double gamma_exact(double & x)
{
    // Split the function domain into three intervals:
    // (0, 0.001), [0.001, 12), and (12, infinity)

    ///////////////////////////////////////////////////////////////////////////
    // First interval: (0, 0.001)
	//
	// For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	// So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	// The relative error over this interval is less than 6e-7.

	const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

    if (x < 0.001)
        return 1.0/(x*(1.0 + gamma*x));

    ///////////////////////////////////////////////////////////////////////////
    // Second interval: [0.001, 12)

	if (x < 12.0)
    {
        // The algorithm directly approximates gamma over (1,2) and uses
        // reduction identities to reduce other arguments to this interval.

		double y = x;
        int n = 0;
        bool arg_was_less_than_one = (y < 1.0);

        // Add or subtract integers as necessary to bring y into (1,2)
        // Will correct for this below
        if (arg_was_less_than_one)
        {
            y += 1.0;
        }
        else
        {
            n = static_cast<int> (floor(y)) - 1;  // will use n later
            y -= n;
        }

        // numerator coefficients for approximation over the interval (1,2)
        static const double p[] =
        {
            -1.71618513886549492533811E+0,
             2.47656508055759199108314E+1,
            -3.79804256470945635097577E+2,
             6.29331155312818442661052E+2,
             8.66966202790413211295064E+2,
            -3.14512729688483675254357E+4,
            -3.61444134186911729807069E+4,
             6.64561438202405440627855E+4
        };

        // denominator coefficients for approximation over the interval (1,2)
        static const double q[] =
        {
            -3.08402300119738975254353E+1,
             3.15350626979604161529144E+2,
            -1.01515636749021914166146E+3,
            -3.10777167157231109440444E+3,
             2.25381184209801510330112E+4,
             4.75584627752788110767815E+3,
            -1.34659959864969306392456E+5,
            -1.15132259675553483497211E+5
        };

        double num = 0.0;
        double den = 1.0;
        int i;

        double z = y - 1;
        for (i = 0; i < 8; i++)
        {
            num = (num + p[i])*z;
            den = den*z + q[i];
        }
        double result = num/den + 1.0;

        // Apply correction if argument was not initially in (1,2)
        if (arg_was_less_than_one)
        {
            // Use identity gamma(z) = gamma(z+1)/z
            // The variable "result" now holds gamma of the original y + 1
            // Thus we use y-1 to get back the orginal y.
            result /= (y-1.0);
        }
        else
        {
            // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
            for (i = 0; i < n; i++)
                result *= y++;
        }

		return result;
    }

    ///////////////////////////////////////////////////////////////////////////
    // Third interval: [12, infinity)

    if (x > 171.624)
    {
		// Correct answer too large to display. Force +infinity.
		double temp = DBL_MAX;
		return temp*2.0;
    }

    return exp(lngamma_exact(x));
}


double lngamma_exact (double & x)
{

    if (x < 12.0)
    {
        return log(fabs(gamma_exact(x)));
    }

	// Abramowitz and Stegun 6.1.41
    // Asymptotic series should be good to at least 11 or 12 figures
    // For error analysis, see Whittiker and Watson
    // A Course in Modern Analysis (1927), page 252

    static const double c[8] =
    {
		 1.0/12.0,
		-1.0/360.0,
		1.0/1260.0,
		-1.0/1680.0,
		1.0/1188.0,
		-691.0/360360.0,
		1.0/156.0,
		-3617.0/122400.0
    };
    double z = 1.0/(x*x);
    double sum = c[7];
    for (int i=6; i >= 0; i--)
    {
        sum *= z;
        sum += c[i];
    }
    double series = sum/x;

    static const double halfLogTwoPi = 0.91893853320467274178032973640562;
    double logGamma = (x - 0.5)*log(x) - x + halfLogTwoPi + series;
	return logGamma;
}

double n_choose_k (int n, double k)
{
    if ( n >= 0 && k == 0)
        return 1;
    else if(n == 0 && k >= 0)
        return 0;
    else
        return n_choose_k(n-1, k) + n_choose_k(n-1,k-1);
}

double incomplete_beta (double a, double b, double x)
{
    double Ix = 0;
    for (int i=a; i<(a+b); i++) {
        Ix += n_choose_k(a+b-1,i)*pow(x,i)*pow(1-x,a+b-1-i);
    }
    return Ix;
}

double sgn (double x)
{

    if (x > 0)
        return 1;
    else if(x < 0)
        return -1;
    else
        return 0;
}


double incomplete_gamma (double a, double x)
{
    const int ITMAX = 100;
    const double EPS = 2.22045e-016;
    int n;
    double sum;
    double gamser;
    double gln = randnumbers::lngamma_exact(a);
    double ap = a;
    double del = sum = 1.0/a;
    for (n=0; n<ITMAX; n++) {
        ++ap;
        del *= x/ap;
        sum += del;
        if(fabs(del)<fabs(sum)*EPS) {
            gamser = sum*exp(-x+log(x)-gln);
        }
    }
    return gamser;
}


} // end: namespace randnumbers
