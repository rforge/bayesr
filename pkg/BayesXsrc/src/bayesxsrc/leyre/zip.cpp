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



#include "zip.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop
#endif

using std::ios;

namespace MCMC
{



// Default constructor

DISTRIBUTION_zip::DISTRIBUTION_zip(void)
: DISTRIBUTION()
{
    family = "Zero Inflated Count Data Distributions";
}


// Constructor without offset
DISTRIBUTION_zip::DISTRIBUTION_zip(const double & a,
                             const double & b,
                             const double & pv,const zipvertopt & vo,
                             const zippropscale & psc,bool hie,
                             MCMCoptions * o,
                             const datamatrix & r,
                             const ST::string & p,const ST::string & ps,
                             const datamatrix & w)

   : DISTRIBUTION(o,r,w,p,ps)
{

create(o, a, b, pv, vo, psc, hie, ps);

}

//Constructor with offset

DISTRIBUTION_zip::DISTRIBUTION_zip(const double & a,
                                    const double & b,
                                    const double & pv,const zipvertopt & vo,
                                    const zippropscale & psc,bool hie,
                                    const datamatrix & offset,
                                    MCMCoptions * o,
                                    const datamatrix & r,
                                    const ST::string & p,const ST::string & ps,
                                    const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w,p,ps)
{

create(o, a, b, pv, vo, psc, hie, ps);

}


void DISTRIBUTION_zip::create(MCMCoptions * o, const double & a,
                                    const double & b, const double & pv,
                                    const zipvertopt & vo, const zippropscale & psc,
                                    bool hie, const ST::string & ps)
{

    if(nrobs>500) oversize = true;
    else oversize = false;

    hierarchical = hie;
    family = "Zero Inflated Count Data Distributions";
    if(ver!=zip)
        scaleexisting = true;
    else
        scaleexisting = false;

// Für die Acceptance-Quote der S-P
    accept=datamatrix(nrobs+3 ,1, 0);

    nu = datamatrix(nrobs, 1, 1);
    unsigned l = ps.length();
    ST::string pathnu = ps.substr(0, l-9) + "nu_sample.raw";

    if(oversize)     // Für die KFZ-Daten!!!!
    {
        nusave = FULLCOND(o, datamatrix(nrobs, 1), "Multiplicative Random Effects",
                nrobs, 1, pathnu);
        nusave.setflags(MCMC::norelchange | MCMC::nooutput | MCMC::nosamples);

        ST::string pathnukfz = ps.substr(0, l-9) + "nu_long_sample.raw";
        nusavekfz = FULLCOND(o, datamatrix(10, 1), "Multiplicative Random Effects reduced!",
                10, 1, pathnukfz);
        nusavekfz.setflags(MCMC::norelchange | MCMC::nooutput);
    }
    else  //  Normale Version!!!!!
    {
        nusave = FULLCOND(o, datamatrix(nrobs, 1), "Multiplicative Random Effects",
                nrobs, 1, pathnu);
        nusave.setflags(MCMC::norelchange | MCMC::nooutput);
    }

    if(hierarchical)
    {
        hierint = datamatrix(1, 1, 0);
        ST::string pathhierint = ps.substr(0, l-9) + "hierarchical_intercept_sample.raw";
        hierintsave = FULLCOND(o, datamatrix(1, 1), "Hierarchical intercept",
                1, 1, pathhierint);
        hierintsave.setflags(MCMC::norelchange | MCMC::nooutput);
    }


    pvar = datamatrix(nrobs+3,1,0.5);
    a_pri = a;
    b_pri = datamatrix(1, 1, b);
    ST::string pathb = ps.substr(0, l-9) + "b_sample.raw";
    b_pri_save = FULLCOND(o, datamatrix(1, 1), "b-Hyperparameter for scale",
            1, 1, pathb);
    b_pri_save.setflags(MCMC::norelchange | MCMC::nooutput);


    theta = datamatrix(1, 1, 0);
    ST::string paththeta = ps.substr(0, l-9) + "theta_sample.raw";
    theta_save = FULLCOND(o, datamatrix(1, 1), "Probability for zero count", 1, 1, paththeta);
    theta_save.setflags(MCMC::norelchange | MCMC::nooutput);

    //Number of zero-observations in the data (Helpvariable)

    m = datamatrix(1, 1, 0);


    prop_var = pv;
    ver = vo;
    pscale = psc;

    sum_nu = datamatrix(1, 1, 0);

    sum2_nu = datamatrix(1, 1, 0);



}

DISTRIBUTION_zip::DISTRIBUTION_zip(const DISTRIBUTION_zip &nd)
                          : DISTRIBUTION(DISTRIBUTION(nd))
{

    oversize = nd.oversize;
    hierarchical = nd.hierarchical;
    accept = nd.accept;
    nu = nd.nu;
    nusave = nd.nusave;
    nusavekfz = nd.nusavekfz;
    hierint = nd.hierint;
    hierintsave = nd.hierintsave;
    pvar = nd.pvar;
    a_pri = nd.a_pri;
    b_pri = nd.b_pri;
    b_pri_save = nd.b_pri_save;
    theta = nd.theta;
    theta_save = nd.theta_save;
    m = nd.m;
    prop_var = nd.prop_var;
    ver = nd.ver;
    pscale = nd.pscale;
    sum_nu = nd.sum_nu;
    sum2_nu = nd.sum2_nu;

}


const DISTRIBUTION_zip &
   DISTRIBUTION_zip::operator=(const DISTRIBUTION_zip & nd)
{

     if (this==&nd)
	   return *this;
    DISTRIBUTION::operator=(DISTRIBUTION(nd));

    oversize = nd.oversize;
    hierarchical = nd.hierarchical;
    accept = nd.accept;
    nu = nd.nu;
    nusave = nd.nusave;
    nusavekfz = nd.nusavekfz;
    hierint = nd.hierint;
    hierintsave = nd.hierintsave;
    pvar = nd.pvar;
    a_pri = nd.a_pri;
    b_pri = nd.b_pri;
    b_pri_save = nd.b_pri_save;
    theta = nd.theta;
    theta_save = nd.theta_save;
    m = nd.m;
    prop_var = nd.prop_var;
    ver = nd.ver;
    pscale = nd.pscale;
    sum_nu = nd.sum_nu;
    sum2_nu = nd.sum2_nu;

    return *this;
}


// Log_Likelihood
double DISTRIBUTION_zip::loglikelihood(double * response,
                       double * linpred,double * weight,const int & i) const
{

    // For the ZINB model.

  if(ver==zinb)
  {
    if(*response==0)
      return log(theta(0,0) + (1-theta(0,0)) * pow(scale(0,0)/(scale(0,0)+exp(*linpred)), scale(0,0)));
    else
      return - (*response + scale(0,0))*log(exp(*linpred) + scale(0,0)) + *response*(*linpred);

  }
  // if ZIP, ZIPGA or ZIPIG

  else
  {
    if(*response==0)
      return log(theta(0,0) + (1-theta(0,0)) * exp(-exp(*linpred)));
    else
      return *weight * (*response * *linpred - exp(*linpred));
  }
}


void DISTRIBUTION_zip::compute_mu(const double * linpred,
  double * mu) const
{

    *mu = exp(*linpred);

}


void DISTRIBUTION_zip::compute_mu_notransform(
const double * linpred,double * mu) const
  {
    *mu = exp(*linpred);
  }


void DISTRIBUTION_zip::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale, const int & i) const
{

    // For ZINB models!
    // Diferenciate between zero and nonzero observation.

    if(ver==zinb)
    {
        double h = scale(0,0)/(scale(0,0)+*mu);

        // Zero response!
        // scale = s, theta = z
        // DEV = -2*log(P_ZINB(y = 0|...)) =
        //     = -2*log(z + (1-z)*(s/(s+mu))^s)

        if(*response==0)
        {
            *deviance= -2*log(theta(0,0) + (1-theta(0,0))*pow(h, scale(0,0)));
            *deviancesat = *deviance;
        }
        // Nonzero response!
        // scale = s, theta = z
        // DEV = -2*log(P_ZINB(y|...)) =
        //     = -2*(log(1-z)+logG(s+y)-logG(y+1)-logG(s)+s*log(s)-s*log(s+mu)
        //        +y*log(mu)-y*log(s+mu))

        else
        {
            *deviance= -2*(log(1-theta(0,0))+lgamma(scale(0,0)+*response)-
                lgamma(scale(0,0))-lgamma(*response+1)+
                scale(0,0)*log(h)+ *response*log(1-h));
            *deviancesat = 2*(*response * log(*response/ *mu)+
                (scale(0,0)+*response)*log((scale(0,0)+*mu)/
                (scale(0,0)+*response)));
        }
    }

    // For ZIP, ZIPGA and ZIPIG models!!
    // Diferenciate between zero and nonzero observation.

    else // zip, zipga, zipig
    {
        if(*response==0)
        {

        // Zero response!
        // theta = z
        // DEV = -2*log(P_ZIP(y = 0|...)) =
        //     = -2*log(z + (1-z)*exp(mu))

            *deviance = -2*log(theta(0,0)+(1-theta(0,0))*exp(-*mu));
            *deviancesat = *deviance;
        }
        else
        {

        // Nonzero response!
        // theta = z
        // DEV = -2*log(P_ZIP(y|...)) =
        //     = -2*(log(1-z) - mu + y*log(mu) - log(y!))

            *deviance = -2*(log(1-theta(0,0))-*mu+*response*log(*mu)-lgamma(*response+1));
            *deviancesat = -2*(*response-*mu+*response*log(*mu/ *response));
        }

    }
}


void  DISTRIBUTION_zip::tilde_y(datamatrix & tildey,datamatrix & m,const unsigned & col,
                                     const bool & current,const datamatrix & w)
 {
  unsigned i;

  double * splinework = m.getV();
  double * weightwork = w.getV();
  double * respwork = response.getV();
  double * linwork = (*linpred_current).getV();
  double * ywork = tildey.getV();
  double * thetawork = theta.getV();
  double * scalework = scale.getV();
  double l, help;
  double explin;

  for (i=0;i<nrobs;i++,splinework++,ywork++,respwork++,weightwork++, linwork++)
  {
    explin = exp(*linwork);
    l = loglikelihood(respwork, linwork, weightwork, i);

    if(ver==zinb)
    {
        if(*respwork==0)
        {
            help = *scalework/(*scalework+ explin);
            *ywork =  -l/(help*(l-*thetawork*explin));
        }
        else // *respwork !=0
            *ywork = (*respwork - explin)/explin;
    } // end zinb
    else // zip zipga zipig
    {
        if(*respwork==0)
        {
            *ywork = -l/(l-*thetawork*explin);
        }
        else // *respwork !=0
        {
            *ywork = (*respwork - explin)/explin;
        }
    } // end zip zipga zipig

      *ywork += *splinework;
   } //end for...

}





double DISTRIBUTION_zip::compute_weight(double * linpred,double *weight,
                        const int & i,const unsigned & col) const
{

    double * respwork = response.getV();
    respwork += i;
    double l = loglikelihood(respwork, linpred, weight, i);
    l = exp(l);
    double explin = exp(*linpred);
    double *scalework = scale.getV();
    double *thetawork = theta.getV();
    double help;


    if(optionsp->get_nriter()<1) // Für Posteriormode!!!!!  ohne Zero Inflated!!!!!!
    {
        if(ver!=zip)
            return *weight * explin* *scalework/(explin + *scalework);
        else
            return *weight * explin;
    }
    else  // Für MCMC!!!!!!
    {
        if(ver==zinb)
        {
            help = *scalework/(*scalework+ explin);
            if(*respwork==0)
            {
                help = pow(help, *scalework+2);
                return *weight*explin*(1-*thetawork)*help*(l-*thetawork*explin)/(l*l);
            }
            else return *weight * explin* help;
        }
        else // if zip, zipga, zipig
        {
            if(*respwork==0) return *weight*explin*exp(-explin)*(1-*thetawork)*(l-*thetawork*explin)/(l*l);
            else return *weight * explin;
        }
    }
}



double DISTRIBUTION_zip::compute_gmu(double * linpred,
                                          const unsigned & col) const
{
    return 1.0/exp(*linpred);
}


double DISTRIBUTION_zip::compute_IWLS(double * response,double * linpred,
                      double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col)
  {


  double l = loglikelihood(response, linpred, weight, i);
  l = exp(l);
  double explin = exp(*linpred);
  double *scalework = scale.getV();
  double *thetawork = theta.getV();
  double help;
  double help2;


        if(ver==zinb)
        {
            if(*response==0)
            {
                help = *scalework/(*scalework+ explin);
                *tildey = -l/(help*(l-*thetawork*explin));
                if (weightyes)
                {
                    help2 = pow(help, *scalework+2);
                    *weightiwls = *weight*explin*(1-*thetawork)*help2*(l-*thetawork*explin)/(l*l);
                }
                return log(*thetawork+(1-*thetawork)*pow(help, *scalework));
            }
            else //  response != 0 !!!
            {
                *tildey = (*response - explin)/explin;
                if (weightyes)
                    *weightiwls = *weight*explin* *scalework/(explin + *scalework);
                return -(*response+*scalework)*log(explin + *scalework) + *response*(*linpred);
            }
        } // end zinb
        else // if zip, zipga, zipig
        {
            if(*response==0)
            {
                *tildey = -l/(l-*thetawork*explin);
                help = exp(-explin);
                if(weightyes)
                    *weightiwls = *weight*explin*(1-*thetawork)*help*(l-*thetawork*explin)/(l*l);
                return log(*thetawork + (1-*thetawork)*help);
            }
            else // response !=0 !!!!
            {
              *tildey =  (*response - explin)/explin;
              if (weightyes)
                    *weightiwls = *weight * explin;
              return *weight * (*response * *linpred - explin);
            }
        } //end zip zipga zipig
  }


void DISTRIBUTION_zip::compute_IWLS_weight_tildey(
                              double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {
  double l = loglikelihood(response, linpred, weight, i);
  l = exp(l);
  double explin = exp(*linpred);
  double *scalework = scale.getV();
  double *thetawork = theta.getV();
  double help;
  double help2;


        if(ver==zinb)
        {
            if(*response==0)
            {
                help = *scalework/(*scalework+ explin);
                *tildey = -l/(help*(l-*thetawork*explin));
                help2 = pow(help, *scalework+2);
                *weightiwls = *weight*explin*(1-*thetawork)*help2*(l-*thetawork*explin)/(l*l);
            }
            else //  response != 0 !!!
            {
                *tildey = (*response - explin)/explin;
                *weightiwls = *weight*explin* *scalework/(explin + *scalework);
            }
        } // end zinb
        else // if zip, zipga, zipig
        {
            if(*response==0)
            {
                *tildey = -l/(l-*thetawork*explin);
                help = exp(-explin);
                *weightiwls = *weight*explin*(1-*thetawork)*help*(l-*thetawork*explin)/(l*l);
            }
            else // response !=0 !!!!
            {
              *tildey =  (*response - explin)/explin;
              *weightiwls = *weight * explin;
            }
        } //end zip zipga zipig
  }



void DISTRIBUTION_zip::outoptions(void)
{
    DISTRIBUTION::outoptions();
    ST::string dist;
    if(ver == zinb) dist = "zero inflated negative binomial";
    else if(ver == zip) dist = "zero inflated Poisson";
    else if(ver == zipga) dist = "zero inflated poisson-gamma";
    else dist = "zero inflated poisson-inverse gaussian";
    ST::string proposal;
    if(pscale == gamzip) proposal = "gamma";
    else proposal = "uniform";
// include program code here
    optionsp->out("  Options for the response variable: \n");
    optionsp->out("\n");
    optionsp->out("     Distribution: " + dist + "\n");
    optionsp->out("\n");
    optionsp->out("  Options for the scale parameter: \n");
    optionsp->out("\n");
    optionsp->out("     Proposal distribution: " + proposal + "\n");
    optionsp->out("     Hyperparameter a for the prior: " +
                    ST::doubletostring(a_pri,6) + "\n");
    optionsp->out("\n");
}

void DISTRIBUTION_zip::outresults(void)
{

    DISTRIBUTION::outresults();
    optionsp->out("\n\n");

    unsigned i, l = pathresultsscale.length();

    double lower1 = Scalesave.get_lower1();
    double lower2 = Scalesave.get_lower2();
    double upper1 = Scalesave.get_upper1();
    double upper2 = Scalesave.get_upper2();


    ST::string l1 = ST::doubletostring(lower1,4);
    ST::string l2 = ST::doubletostring(lower2,4);
    ST::string u1 = ST::doubletostring(upper1,4);
    ST::string u2 = ST::doubletostring(upper2,4);

    ST::string nl1 = l1.replaceallsigns('.','p');
    ST::string nl2 = l2.replaceallsigns('.','p');
    ST::string nu1 = u1.replaceallsigns('.','p');
    ST::string nu2 = u2.replaceallsigns('.','p');


    if(ver!=zip)
    {
        double sum=0.0;

        if(ver!=zinb)
        {
            if(ver==zipig)
            {
                double *acceptwork=accept.getV();
                unsigned i;
                acceptwork +=1;
                for(i=0; i<(accept.rows()-2); i++, acceptwork++)
                {
                    sum += *acceptwork;
                }
                sum = sum/nrobs;
                sum = optionsp->get_nriter()-optionsp->get_burnin();

                optionsp->out("\n\n");
                optionsp->out("  Acceptance rate for the parameter block nu:   " +
                ST::doubletostring(double(sum)/
                    double(optionsp->get_nriter()-optionsp->get_burnin())*100, 4)
                    + " %" + "\n");
                optionsp->out("\n");
            }
            else
            {
                optionsp->out("\n\n");
                optionsp->out("  Acceptance rate for the parameter block nu:   " +
                ST::doubletostring(100) + " %" + "\n");
                optionsp->out("\n");
            }
    }





//Outfile for the nu-sample
//Write out nu



    if(ver != zinb)
    {
        if(!oversize) // Für die Normale Version!!!!
        {
            nusave.outresults();
            ST::string pathnu = pathresultsscale.substr(0, l-10) + "_nu_sample.raw";
            nusave.get_samples(pathnu);
        }
        else    // Für die KFZ-Daten!!!!!
        {
            nusavekfz.outresults();
            ST::string pathnukfz = pathresultsscale.substr(0, l-10) + "_nu_sample.raw";
            nusavekfz.get_samples(pathnukfz);
        }

        ST::string pathnures = pathresultsscale.substr(0, l-10) + "_nu.res";


        if(!oversize) // für die normale Version!!!!!
        {
            double * workmean = nusave.get_betameanp();
            double * workstddev = nusave.get_betavarp();
            double * workbetaqu50 = nusave.get_betaqu50p();
            double * workbetaqu_l1_lower_p = nusave.get_beta_lower1_p();
            double * workbetaqu_l2_lower_p = nusave.get_beta_lower2_p();
            double * workbetaqu_l1_upper_p = nusave.get_beta_upper1_p();
            double * workbetaqu_l2_upper_p = nusave.get_beta_upper2_p();

            ofstream outnu(pathnures.strtochar());
            outnu << "nu" << "   ";
            outnu << "pmean" << "   ";
            outnu << "stddev" << "   ";
            outnu << "pqu"  << l1  << "   ";
            outnu << "pqu"  << l2  << "   ";
            outnu << "pmed   ";
            outnu << "pqu"  << u1  << "   ";
            outnu << "pqu"  << u2  << "   ";
            outnu << endl;

            for(i=0;i<nrobs;i++,workmean++, workstddev++, workbetaqu_l1_lower_p++,
                           workbetaqu_l2_lower_p++,workbetaqu50++,
                           workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
            {
                outnu << (i+1) << "   ";
                outnu << *workmean << "   ";
                outnu << *workstddev << "   ";
                outnu << *workbetaqu_l1_lower_p << "   ";
                outnu << *workbetaqu_l2_lower_p << "   ";
                outnu << *workbetaqu50 << "   ";
                outnu << *workbetaqu_l2_upper_p << "   ";
                outnu << *workbetaqu_l1_upper_p << "   ";
                outnu << endl;

            }
        }
        else   // Für KFZ-Daten!!!
        {
            double * workmean = nusave.get_betameanp();
            double * workstddev = nusave.get_betavarp();

            ofstream outnu(pathnures.strtochar());
            outnu << "nu" << "   ";
            outnu << "pmean" << "   ";
            outnu << "stddev" << "   ";
            outnu << endl;

            for(i=0;i<nrobs;i++,workmean++,workstddev++)
            {
            outnu << (i+1) << "   ";
            outnu << *workmean << "   ";
            outnu << *workstddev << "   ";
            outnu << endl;
            }
        }

        if(hierarchical)
        {
            hierintsave.outresults();
            ST::string pathhier = pathresultsscale.substr(0, l-10) + "_hierarchical_intercept_sample.raw";
            hierintsave.get_samples(pathhier);

            ST::string pathhierintres = pathresultsscale.substr(0, l-10) + "_hierarchical_intercept.res";

            ofstream outhierint(pathhierintres.strtochar());

            outhierint << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                    "   pmed   pqu" << nu1 << "   pqu" << nu2 << endl;

            outhierint << hierintsave.get_betamean(0,0) << "  " <<
                    sqrt(hierintsave.get_betavar(0,0)) << "  "<<
                    hierintsave.get_beta_lower1(0,0) << "  "<<
                    hierintsave.get_beta_lower2(0,0) << "  "<<
                    hierintsave.get_betaqu50(0,0) << "  "<<
                    hierintsave.get_beta_upper2(0,0) << "  "<<
                    hierintsave.get_beta_upper1(0,0) << endl;
        }

    }

//Outfile for the b-pri-sample
//Write out b-pri

    b_pri_save.outresults();
    ST::string pathb = pathresultsscale.substr(0, l-10) + "_b_sample.raw";
    b_pri_save.get_samples(pathb);


    ST::string pathbres = pathresultsscale.substr(0, l-10) + "_b_pri.res";

    ofstream outb(pathbres.strtochar());
    outb << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                    "   pmed   pqu" << nu1 << "   pqu" << nu2 << endl;

    outb << b_pri_save.get_betamean(0,0) << "  " <<
                    sqrt(b_pri_save.get_betavar(0,0)) << "  "<<
                    b_pri_save.get_beta_lower1(0,0) << "  "<<
                    b_pri_save.get_beta_lower2(0,0) << "  "<<
                    b_pri_save.get_betaqu50(0,0) << "  "<<
                    b_pri_save.get_beta_upper2(0,0) << "  "<<
                    b_pri_save.get_beta_upper1(0,0) << endl;


//Outfile for the theta-sample
//Write out theta

    } // end if ver!=zip

    optionsp->out("\n\n");
    optionsp->out("  Acceptance rate for the parameter theta:   " +
    ST::doubletostring(double(accept(nrobs+2, 0))/
    double(optionsp->get_nriter()-optionsp->get_burnin())*100, 4) + " %" + "\n");
    optionsp->out("\n");


    theta_save.outresults();
    ST::string paththeta = pathresultsscale.substr(0, l-10) + "_theta_sample.raw";
    theta_save.get_samples(paththeta);


    ST::string paththetares = pathresultsscale.substr(0, l-10) + "_theta.res";

    ofstream outtheta(paththetares.strtochar());
    outtheta << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                    "   pmed   pqu" << nu1 << "   pqu" << nu2 << endl;

    outtheta << theta_save.get_betamean(0,0) << "  " <<
                    sqrt(theta_save.get_betavar(0,0)) << "  "<<
                    theta_save.get_beta_lower1(0,0) << "  "<<
                    theta_save.get_beta_lower2(0,0) << "  "<<
                    theta_save.get_betaqu50(0,0) << "  "<<
                    theta_save.get_beta_upper2(0,0) << "  "<<
                    theta_save.get_beta_upper1(0,0) << endl;





    if(ver!=zinb && ver!=zip)
    {
        unsigned j;
        ST::string hilfe;
        double *pwork = pvar.getV();
        double *acceptwork=accept.getV();
        pwork+=1;
        acceptwork+=1;
        ofstream out;
        hilfe = "_pvar.raw";
        unsigned l = pathresultsscale.length();
        ST::string pathp = pathresultsscale.substr(0, l-10) + hilfe;
        out.open(pathp.strtochar(), ios::out);
        out << "pvar" << " " << "accept" << endl;

        for(j=1; j < pvar.rows(); j++, pwork++, acceptwork++)
        {
            out << j << " " << *pwork << " " << *acceptwork << endl;
        }
    }

}





void DISTRIBUTION_zip::update(void)
{

    // For hierarchical models, initialize the multiplicative random effects in
    // the first iteration with exp(intercept). And for all the models initialize
    // the zero inflation parameter theta with the value of proportion of zero
    // values in the data. Is a good approximation only for Distributions with
    // a rather large mean...

    if(optionsp->get_nriter()==1)
    {
     unsigned i;
     double *nuwork = nu.getV();
	 double *respwork = response.getV();
	 double *mwork = m.getV();

     // Initialize Hierarchical intercept

	 if(hierarchical)
	 {
	   double *intwork = hierint.getV();
	   double expint = exp(*intwork);
       for(i=0; i<nrobs; i++, nuwork++)
       {
            *nuwork = expint;
       }
	 }

     // Initialize theta

     for(i=0; i<nrobs; i++, respwork++)
     {
	   if(*respwork == 0) *mwork += 1;
     }
     theta(0,0) = *mwork/nrobs;
     *mwork = nrobs - *mwork;
    }

    // Calling update functions!


    // First only models with overdispersion

    if(ver!=zip)
    {

        // Models with multiplicative random effects
        // Update of the multiplicative random effects
        // If hierarchical = true, then update of hierarchical intercept.

        if(ver!=zinb)       //PoGa oder PoIg
        {
            update_nu();
            nusave.update();
            if(oversize) nusavekfz.update();
            if(hierarchical)
            {
                update_hierint();
                hierintsave.update();
            }
        }

        // Update scale for all the models with overdispersion. Differentiations
        // for the models are included in the function.


        update_scale();

        // Update of the hyperparameter b of the prior for the scale parameter.
        // It is the same function for all the models.

        b_pri(0,0) = update_b_pri();
        b_pri_save.update();
    }

    //  Update theta for all the models!

    update_theta();

    theta_save.update();

    // Calculate the acceptance rates with the help of the vector acceptance
    // when the running is ended.

    if(optionsp->get_nriter()==optionsp->get_iterations())
    {

        if(ver!=zipig)
        {
            acceptancescale = double(accept(0,0))/
                double(optionsp->get_nriter()-optionsp->get_burnin())*100;
        }
        else
        {
            acceptancescale=100;
        }
    }



    DISTRIBUTION::update();
}


bool DISTRIBUTION_zip::posteriormode(void)
  {


   if(!oversize)
    {
        double scaleold = scale(0,0);
        double * respwork = response.getV();
        double * linwork = (*linpred_current).getV();
        double sum = 0;
        double h;
        unsigned i;

        for(i=0; i<nrobs; i++, respwork++, linwork++)
        {
            h = exp(*linwork);
            sum += (*respwork)*(*respwork)/(h*h)-(2*(*respwork)+1)/h;
        }

        scale(0,0) = nrobs/(sum+nrobs);

        double norm = scaleold-scale(0,0);
        norm = (norm*norm)/scaleold*scaleold;

        if (norm <= 0.00001)
            return true;
        else
            return false;
     }
     else
        return true;
}


bool DISTRIBUTION_zip::posteriormode_converged(const unsigned & itnr)
{
   if(!oversize)
    {
        if (itnr > 1)
        {
            if ( fabs(scale(0,0)-scale_mode(0,0))/scale_mode(0,0) < 0.00001)
              return true;
            else
              return false;
        }
        return false;
    }
    else
        return true;
}



double DISTRIBUTION_zip::update_nu(void)
{
    unsigned i;
    double *nuwork = nu.getV();
    double *respwork = response.getV();
    double *linwork = (*linpred_current).getV();
    double *scalework=scale.getV();
    double mu;
    double nuaux;   // auxiliary nu
    double * nusavep = nusave.getbetapointer();
    double *sum=sum_nu.getV();
    *sum=0.0;
    double *sum2=sum2_nu.getV();
    *sum2=0.0;
    double log_proposal_ratio;
    double log_likelihood_ratio;
    double log_priori_ratio;
    double *thetawork = theta.getV();
    double alpha;
    double h;
    double *acceptwork=accept.getV();
    acceptwork++;

//    if(hierarchical)
    double *intwork = hierint.getV(); // only for hierarchical
    double expint = exp(*intwork);   // only for hierarchical


    //  When oversize=True, then only few samplings for the multiplicative
    // random effects are stored.

    double * nusavekfzp = nusavekfz.getbetapointer();
    unsigned helpsum = (nrobs-nrobs%10)/10;
    unsigned help;
        if(helpsum%2 == 0) help=helpsum/2;
        else help = (helpsum-1)/2;
    unsigned j = 0;



    //Update nu: M-H algorithm for all the models!! the only difference in the
    // update step is given by the priori ratio of the corresponding model
    // selection. We differentiate in the implementation between response = 0
    // and response != 0.



	for (i=0;i<nrobs;i++,nuwork++, respwork++,linwork++, nusavep++, acceptwork++)
	{

		mu = exp(*linwork)/(*nuwork);
		nuaux = *nuwork;
		log_proposal_ratio = proposal_nu(i); // nuaux=old value, nuwork=new value

        // Zero response!!

		if(*respwork==0)
		{
			log_likelihood_ratio = log((*thetawork+(1-*thetawork)*exp(-*nuwork*mu))/(*thetawork+(1-*thetawork)*exp(-exp(*linwork))));
		}

        // Nonzero response!!

		else // response != 0
		{
		    	log_likelihood_ratio = mu * (nuaux - *nuwork) + *respwork* (log(*nuwork) - log(nuaux));
		}

        // ZIPGA model = priori of the nu is a Gamma(scale, scale)
        // and a Gamma(scale, scale/mu) for the hierarchical version.

		if(ver==zipga)
		{
			if(hierarchical) log_priori_ratio = (*scalework-1)*log(*nuwork/nuaux)+*scalework*(nuaux-*nuwork)/expint;
			else log_priori_ratio = (*scalework-1)*log(*nuwork/nuaux)+*scalework*(nuaux-*nuwork);
		}

        // ZIPIG model = priori of the nu is a IGauss(1, scale)
        // and a IGauss(mu, mu*scale) for the hierarchical version.

		else // poig
		{
		    	if(hierarchical) log_priori_ratio = -1.5*log(*nuwork/nuaux)+*scalework*0.5*((nuaux-*nuwork)/expint+expint*(1/nuaux-1/(*nuwork)));
		    	else log_priori_ratio = -1.5*log(*nuwork/nuaux)+*scalework*0.5*(nuaux +1/nuaux-*nuwork - 1/(*nuwork));
		}

		alpha = log_likelihood_ratio + log_priori_ratio + log_proposal_ratio;
	        h =log(randnumbers::uniform( ));
	        if(h<=alpha)
	        {
		   *linwork += log(*nuwork/(nuaux));
	           *acceptwork=*acceptwork+1;
	        }
	        else   *nuwork = nuaux;
		if(ver==zipga)
		{
		   *sum += *nuwork;
		   *sum2 += log(*nuwork);
		}
		else // poig...
		{
	           *sum += *nuwork;
	           *sum2 += 1/(*nuwork);
		}
	        *nusavep = *nuwork;
	        if(oversize)            // For KFZ_Data!!!!
	        {
	                if(i==help && j<10)
                	{
	                    *nusavekfzp=*nuwork;
	                    help += helpsum;
	                    nusavekfzp++;
	                    j++;
	                }
	        }

            // Tuning function for the acceptance rates of the M-H algorithm
            // It will be called only in the burnin phase and each 100 iterations.

            if(optionsp->get_nriter() % 100 == 0 &&
	            optionsp->get_nriter() <= optionsp->get_burnin())
            {
                pwork_tunin(i+1);
            }
	}

    return 0;

}


double DISTRIBUTION_zip::update_hierint(void)
{
    double *intwork = hierint.getV();
    double * hierintsavep = hierintsave.getbetapointer();

    *hierintsavep = *intwork;

    // Tuning function for the acceptance rates of the M-H algorithm
    // It will be called only in the burnin phase and each 100 iterations.

    if(optionsp->get_nriter() % 100 == 0 &&
            optionsp->get_nriter() <= optionsp->get_burnin())
    {
        pwork_tunin(nrobs + 1);
    }
    return 0;
}


double DISTRIBUTION_zip::update_scale(void) const
{
    double *scalework = scale.getV();

    // The update for the ZINB and ZIPGA models is a M-H step!

    if(ver!=zipig)
    {
        double scaleaux=*scalework;
        double log_likelihood_ratio;
        double priori_ratio;
        double proposal_ratio;
        double alpha;
        double *acceptwork=accept.getV();

        // Proposal value for scale depending on the proposal distribution
        // And proposal_ratio also depending on the proposal distribution

        proposal_ratio = proposal_scale();

        // After this, *scalework stores the new value and scaleaux the old one

        // For a ZINB model:

        if(ver==zinb)  log_likelihood_ratio = log_nbin(*scalework, scaleaux);

        // For a ZIPGA model:

        else    // ver==poga
        {
         if(hierarchical) log_likelihood_ratio = log_gamma_likelihood_hier(scaleaux, *scalework);
         else log_likelihood_ratio = log_gamma_likelihood(scaleaux, *scalework);
        }

        priori_ratio = (a_pri-1)*(log(*scalework)-log(scaleaux))+
                         b_pri(0,0)*(scaleaux - *scalework);

        alpha=priori_ratio + log_likelihood_ratio + proposal_ratio;
        double h =log(randnumbers::uniform( ));
        if(h<=alpha) *acceptwork=*acceptwork+1;
        else *scalework = scaleaux;

        // Tuning function for the acceptance rates of the M-H algorithm
        // It will be called only in the burnin phase and each 100 iterations.

        if(optionsp->get_nriter() % 100 == 0 &&
            optionsp->get_nriter() <= optionsp->get_burnin())
        {
            pwork_tunin(0);
        }
    }

    // The update for the POIG model is a Gibbs-step!
    // Formula (5.21) from diss

    else    // if ver==poig
    {
        double *sum = sum_nu.getV();
        double *sum2 = sum2_nu.getV();

        if(hierarchical)
        {
            double *intwork = hierint.getV();
            double expint = exp(*intwork);
            *scalework = randnumbers::rand_gamma(a_pri+nrobs/2,
                        b_pri(0,0) +0.5*(*sum/expint + *sum2*expint)-nrobs);

        }
        else
        {

            *scalework = randnumbers::rand_gamma(a_pri+nrobs/2,
                        b_pri(0,0) +0.5*(*sum + *sum2)-nrobs);
        }
    }

    return 0;
}


double DISTRIBUTION_zip::update_theta(void)
{

  double *thetawork = theta.getV();
  double * thetasavep = theta_save.getbetapointer();
  double thetaaux = *thetawork;        // For the old value
  double log_likelihood_ratio;
  double log_priori_ratio;
  double log_proposal_ratio;
  double *acceptwork=accept.getV();
  acceptwork += nrobs + 2;

  // Proposal value for theta and proposal_ratio. See below.

  log_proposal_ratio = proposal_theta();

  // After this function, *thetawork stores the proposed value for theta!

  log_priori_ratio = 0.0;

  // Differentiate for the likelihood: ZINB and rest.

  if(ver==zinb)  log_likelihood_ratio=likelihood_zinb(thetaaux);
  else log_likelihood_ratio=likelihood_zirest(thetaaux);

  double alpha=log_priori_ratio + log_likelihood_ratio + log_proposal_ratio;
  double h =log(randnumbers::uniform( ));
  if(h<=alpha)
    *acceptwork=*acceptwork+1;
  else
    *thetawork = thetaaux;

  *thetasavep = *thetawork;

  // Tuning function for the acceptance rates of the M-H algorithm
  // It will be called only in the burnin phase and each 100 iterations.

  if(optionsp->get_nriter() % 100 == 0 &&
            optionsp->get_nriter() <= optionsp->get_burnin())
  {
    pwork_tunin(nrobs+2);
  }


 return 0;

}


double DISTRIBUTION_zip::update_b_pri(void)
{
    double alpha = 1;
    double beta = 0.005;
    double * scalework = scale.getV();
    double * bwork = b_pri.getV();
    double * bp = b_pri_save.getbetapointer();

    // Update by Gibbs-step!
    // Formula (5.22) from diss

    *bwork = randnumbers::rand_gamma(alpha + a_pri, *scalework + beta);
    *bp = *bwork;

    return *bwork;
}


double DISTRIBUTION_zip::proposal_scale(void) const
{


    double *scalework = scale.getV();
    double scaleaux = *scalework;     //alter Wert!!!!!
    double proposal_ratio;
    double *pwork = pvar.getV();


        // For a uniform proposal: pwork stores the width of the uniform proposal.
        // The central point is given by the old value of scale.
        // Formula (5.38) from diss.

       if(pscale == unifzip)       // if uniform
        {
            if(*scalework > *pwork)
            {
                *scalework = scaleaux - *pwork +
                    2*(*pwork)*(randnumbers::uniform());
                if(*scalework > *pwork) proposal_ratio = 1.0;
                else proposal_ratio = (2*(*pwork))/(*scalework + *pwork);
            }
            else
            {
                *scalework = (scaleaux + *pwork)*(randnumbers::uniform());
                if(*scalework < *pwork) proposal_ratio = (scaleaux + *pwork)/
                                                        (*scalework + *pwork);
                else proposal_ratio = (scaleaux + *pwork)/(2*(*pwork));
            }
        proposal_ratio= log(proposal_ratio);
        }

        // For a gamma proposal: pwork stores the variance of the gamma proposal.
        // The mean is given by the old value of scale.
        // Formula (5.39) from diss.

        else     // if gamma
        {
            double a, b, a_neu;//, b_neu;
            a = *scalework * (*scalework)/(*pwork);
            b = *scalework/(*pwork);
            *scalework=randnumbers::rand_gamma(a, b);
            a_neu = *scalework*(*scalework)/(*pwork);
            while(a_neu < exp(-16*log(10.0)))
            {
                *scalework=randnumbers::rand_gamma(a, b);
                a_neu = *scalework*(*scalework)/(*pwork);
            }

            double logsca = log(*scalework);
            double logscaaux = log(scaleaux);
            proposal_ratio = (a_neu-a)*(logscaaux+logsca-log(*pwork))
                            - logscaaux + logsca + lgamma(a) - lgamma(a_neu);



        }

        return proposal_ratio;


}

double DISTRIBUTION_zip::proposal_nu(unsigned i) const
{
    double *nuwork = nu.getV();
    nuwork += i;
    double nuaux = *nuwork;     //alter Wert!!!!!
    double proposal_ratio;
    double *pwork = pvar.getV();
    pwork += (1+i);

    // For a uniform proposal: pwork stores the width of the uniform proposal.
    // The central point is given by the old value of nu_i.
    // Formula (5.38) from diss.

    if(*nuwork > *pwork)
    {
	*nuwork = nuaux -*pwork + 2*(*pwork)*(randnumbers::uniform());
        if(*nuwork > *pwork) proposal_ratio = 1.0;
        else proposal_ratio = (2*(*pwork))/(*nuwork + *pwork);
    }
    else
    {
        *nuwork = (nuaux + *pwork)*(randnumbers::uniform());
        if(*nuwork < *pwork) proposal_ratio = (nuaux + *pwork)/(*nuwork + *pwork);
        else proposal_ratio = (nuaux + *pwork)/(2*(*pwork));
     }

     return log(proposal_ratio);
}

double DISTRIBUTION_zip::proposal_theta(void) const
{
    double *thetawork = theta.getV();
    double thetaaux = *thetawork;     //alter Wert!!!!!
    double log_proposal_ratio;
    double *pwork = pvar.getV();
    pwork += nrobs+2;

    double a = std::max(0.0, thetaaux-*pwork);
    double b = std::min(thetaaux+*pwork, 1.0);
    double a_prop;
    double b_prop;

    // For a uniform proposal: pwork stores the width of the uniform proposal.
    // The central point is given by the old value of theta.
    // Formula (5.43) from diss.

    *thetawork = a +(b-a)*randnumbers::uniform();
    a_prop = std::max(0.0, *thetawork-*pwork);
    b_prop = std::min(*thetawork+*pwork, 1.0);
    return log_proposal_ratio = log(b-a)-log(b_prop-a_prop);
}



double DISTRIBUTION_zip::pwork_tunin(unsigned i) const
{

    double *pwork = pvar.getV();
    pwork += i;
    double *acceptwork = accept.getV();
    acceptwork += i;
    double hilfe = *acceptwork / 100;

//    if(hilfe < 0.2 && *pwork > 0.0001) *pwork = *pwork * 0.1;
    if(hilfe < 0.2) *pwork = *pwork * 0.1;
    if(0.2 < hilfe && hilfe < 0.3 && *pwork > 0.0001) *pwork = *pwork * 0.5;
    if(0.3 < hilfe && hilfe < 0.4 && *pwork > 0.0001) *pwork = *pwork * 0.8;
    if(0.6 < hilfe && hilfe < 0.7) *pwork = *pwork * 1.3;
    if(hilfe > 0.7 && hilfe < 0.8) *pwork = *pwork * 5;
    if(hilfe > 0.8 && *pwork < 10) *pwork = *pwork * 10;

    *acceptwork = 0;


    return 0;

}



double DISTRIBUTION_zip::log_gamma_likelihood(double &s,
                                                    double &s_neu) const
{
    double summe = 0.0;
    double *sum = sum_nu.getV();
    double *sum2 = sum2_nu.getV();
    summe = nrobs*(s_neu*log(s_neu)-s*log(s)+lgamma(s)-lgamma(s_neu)) +
          (s - s_neu)*(*sum - *sum2);
    return summe;
}

double DISTRIBUTION_zip::log_gamma_likelihood_hier(double &s,
                                                    double &s_neu) const
{
    double summe = 0.0;
    double *sum = sum_nu.getV();
    double *sum2 = sum2_nu.getV();
    double * hierintwork = hierint.getV();

    summe = nrobs*(s_neu*log(s_neu)-s*log(s)+(s-s_neu)*(*hierintwork)+
            lgamma(s)-lgamma(s_neu)) +
            (s_neu - s)*(*sum2 - *sum/exp(*hierintwork));
    return summe;
}


//Loggamma-Funktion

double DISTRIBUTION_zip::lgamma(const double & xx) const
{
    if( xx==1 || xx==2 )
    {
        return 0;
    }
    else
    {
        const double stp = 2.50662827465;
        double x, tmp, ser;

        x = xx-1.0;
        tmp = x+5.5;
        tmp = (x+0.5)*log(tmp)-tmp;
        ser = 1.0 + 76.18009173/(x+1.0) - 86.50532033/(x+2.0) +
        24.01409822/(x+3.0) - 1.231739516/(x+4.0) +
        0.120858003e-2/(x+5.0) - 0.536382e-5/(x+6.0);

        return tmp+log(stp*ser);
    }
}



//Loglikelihood-Differenz einer ZINB:
//log(ZINB(lambda, s_neu, theta)/ZINB(lambda, s, theta))

double DISTRIBUTION_zip::log_nbin(const double & s_neu,
                                        const double & s) const
{
    unsigned i;
    double *respwork = response.getV();
    double *linwork = (*linpred_current).getV();
    double *tw = theta.getV();
    double *mw = m.getV();
    double sum = 0.0;

    for(i=0; i<nrobs; i++, respwork++, linwork++)
    {
      if(*respwork==0) sum += log((*tw + (1-*tw)*pow(s_neu/(s_neu+exp(*linwork)), s_neu))/(*tw + (1-*tw)*pow(s/(s+exp(*linwork)), s)));
      else sum += (*respwork+s)*log(exp(*linwork)+s)-(*respwork+s_neu)*log(exp(*linwork)+s_neu);
    }
        sum += *mw*(lgamma(s)-lgamma(s_neu)+s_neu*log(s_neu)-s*log(s));
        return sum;
}


void DISTRIBUTION_zip::add_nu(double m) const
{
    double *nuwork = nu.getV();
    unsigned i;
    for(i=0; i<nrobs; i++, nuwork++)
    {
     *nuwork *= m;
    }
}

//Loglikelihood-Differenz einer ZINB:
//log(ZINB(lambda, s, t)/ZINB(lambda, s, t_neu))


double DISTRIBUTION_zip::likelihood_zinb(const double & t) const
{

    // t: stores the old value!

    unsigned i;
    double *respwork = response.getV();
    double *scalework = scale.getV();
    double *linwork = (*linpred_current).getV();
    double *thetawork = theta.getV(); // stores the new value!!
    double *mwork = m.getV();
    double sum = 0.0;
    double s;
    double h1;


    for(i=0; i<nrobs; i++, respwork++, linwork++)
    {
        if(*respwork ==0)
    	{
              s = *scalework/(*scalework + exp(*linwork));
        	  h1 = pow(s, *scalework);
              sum += log((*thetawork+(1-*thetawork)*h1)/(t+ (1-t)*h1));
    	}
    }
    sum += *mwork *log((1-*thetawork)/(1-t));
    return sum;
}

  //Loglikelihood-Differenz einer ZIP:
  //log(ZIP(lambda, t)/ZIP(lambda, t_neu))

double DISTRIBUTION_zip::likelihood_zirest(const double & t) const
{

    // t: stores the old value!

    unsigned i;
    double *respwork = response.getV();
    double *linwork = (*linpred_current).getV();
    double *thetawork = theta.getV();    // stores the new value!!
    double *mwork = m.getV();
    double sum = 0.0;
    double h1;


    for(i=0; i<nrobs; i++, respwork++, linwork++)
    {
      if(*respwork ==0)
	{
	  h1 = exp(-exp(*linwork));
	  sum += log((*thetawork+h1*(1-*thetawork))/(t+ h1*(1-t)));
	}
    }
    return sum + *mwork *log((1-*thetawork)/(1-t));
}




} // end: namespace MCMC










