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



#include "nbinomial.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop
#endif

using std::ios;

namespace MCMC
{




// Default constructor

DISTRIBUTION_nbinomial::DISTRIBUTION_nbinomial(void)
: DISTRIBUTION()
{
    family = "negative binomial";
}


// Constructor without offset
DISTRIBUTION_nbinomial::DISTRIBUTION_nbinomial(const double & a,
                             const double & b,
                             const double & pv,const vertopt & vo,
                             const propscale & psc,bool hie,
                             MCMCoptions * o,
                             const datamatrix & r,
                             const ST::string & p,const ST::string & ps,
                             const datamatrix & w)

   : DISTRIBUTION(o,r,w,p,ps)
{

create(o, a, b, pv, vo, psc, hie, ps);

}

//Constructor with offset

DISTRIBUTION_nbinomial::DISTRIBUTION_nbinomial(const double & a,
                                    const double & b,
                                    const double & pv,const vertopt & vo,
                                    const propscale & psc,bool hie,
                                    const datamatrix & offset,
                                    MCMCoptions * o,
                                    const datamatrix & r,
                                    const ST::string & p,const ST::string & ps,
                                    const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w,p,ps)
{

create(o, a, b, pv, vo, psc, hie, ps);

}


void DISTRIBUTION_nbinomial::create(MCMCoptions * o, const double & a,
                                    const double & b, const double & pv,
                                    const vertopt & vo, const propscale & psc,
                                    bool hie, const ST::string & ps)
{

    if(nrobs>500) oversize = true;
    else oversize = false;

    hierarchical = hie;
    family = "negative binomial";
    scaleexisting = true;
// Für die Acceptance-Quote der S-P
    accept=datamatrix(nrobs+2 ,1, 0);

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


    pvar = datamatrix(nrobs+2,1,0.5);
    a_pri = a;
    b_pri = datamatrix(1, 1, b);
    ST::string pathb = ps.substr(0, l-9) + "b_sample.raw";
    b_pri_save = FULLCOND(o, datamatrix(1, 1), "b-Hyperparameter for scale",
            1, 1, pathb);
    b_pri_save.setflags(MCMC::norelchange | MCMC::nooutput);

    prop_var = pv;
    ver = vo;
    pscale = psc;

    sum_nu = datamatrix(1, 1, 0);

    sum2_nu = datamatrix(1, 1, 0);



}

DISTRIBUTION_nbinomial::DISTRIBUTION_nbinomial(const DISTRIBUTION_nbinomial &nd)
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
    prop_var = nd.prop_var;
    ver = nd.ver;
    pscale = nd.pscale;
    sum_nu = nd.sum_nu;
    sum2_nu = nd.sum2_nu;

}


const DISTRIBUTION_nbinomial &
   DISTRIBUTION_nbinomial::operator=(const DISTRIBUTION_nbinomial & nd)
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
    prop_var = nd.prop_var;
    ver = nd.ver;
    pscale = nd.pscale;
    sum_nu = nd.sum_nu;
    sum2_nu = nd.sum2_nu;

    return *this;
}


// Log_Likelihood
double DISTRIBUTION_nbinomial::loglikelihood(double * response,
                       double * linpred,double * weight,const int & i) const
{

    // For the negative binomial model
    if (ver==nb)
    {
        return - (*response + scale(0,0))*log(exp(*linpred) + scale(0,0)) +
                    *response*(*linpred);
    }

    // For the POGA and POIG models: corresponds to a Poisson likelihood.
    else
    {
        return *weight * (*response * *linpred - exp(*linpred));
    }

}


void DISTRIBUTION_nbinomial::compute_mu(const double * linpred,
  double * mu) const
{

    *mu = exp(*linpred);

}


void DISTRIBUTION_nbinomial::compute_mu_notransform(
const double * linpred,double * mu) const
  {
    *mu = exp(*linpred);
  }


void DISTRIBUTION_nbinomial::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale, const int & i) const
{

    // Returns only the deviance of a Negative binomial model.
    // Formula: logarithm of the density of a negative binomial distribution
    // NB(mu, d) = G(y+d)/(G(y+1)G(d)) (d/(d+mu))^d (mu/(d+mu))^y
    // DEV = logG(y+d)-logG(y+1)-logG(d) + d*(log(d)-log(mu+d)) + y*(log(mu)-log(mu+d))
    // DEV_SAT = NB(mu, d) - NB(y, d)

    if(*response==0)
    {
        *deviance= -2*scale(0,0)*log(scale(0,0)/(scale(0,0)+*mu));

        *deviancesat = *deviance;
    }
    else
    {
        *deviance= -2*(lgamma(scale(0,0)+*response)-
                lgamma(scale(0,0))-lgamma(*response+1)+
                scale(0,0)*(log(scale(0,0))-log(scale(0,0)+*mu))+
                *response*(log(*mu)-log(scale(0,0)+*mu)));

        *deviancesat = 2*(*response * log(*response/ *mu)+
                (scale(0,0)+*response)*log((scale(0,0)+*mu)/
                (scale(0,0)+*response)));
    }
}



double DISTRIBUTION_nbinomial::compute_weight(double * linpred,double *weight,
                        const int & i,const unsigned & col) const
{

// Für Posteriormode!!!!!
    if(optionsp->get_nriter()<1)
    {
        return *weight * exp(*linpred)*scale(0,0)/(exp(*linpred) +
                        scale(0,0));
    }

// Für MCMC iterations!!!!!!
    else
    {
        // Negative binomial

        if(ver==nb) return *weight * exp(*linpred)*scale(0,0)/(exp(*linpred) +
                        scale(0,0));
        // POGA or POIG models

        else return *weight * exp(*linpred);
    }
}



double DISTRIBUTION_nbinomial::compute_gmu(double * linpred,
                                          const unsigned & col) const
{
    return 1.0/exp(*linpred);
}


double DISTRIBUTION_nbinomial::compute_IWLS(double * response,double * linpred,
                      double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col)
  {


  double el = exp(*linpred);

  // Negative Binomial model:

  if (ver==nb)
    {

    if (weightyes)
      *weightiwls = *weight * el* scale(0,0)/(el + scale(0,0));

    *tildey = (*response - el)/el;

    return  - (*response + scale(0,0))*log(el + scale(0,0)) +
                *response*(*linpred);

    }

    // POGA and POIG models:

  else
    {

    if (weightyes)
     *weightiwls = *weight * el;

    *tildey =  (*response - el)/el;

    return *weight * (*response * *linpred - el);

    }

  }


void DISTRIBUTION_nbinomial::compute_IWLS_weight_tildey(
                              double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {

  double el = exp(*linpred);

  if (ver==nb)
    {

    *weightiwls = *weight * el*scale(0,0)/(el + scale(0,0));

    }
  else
    {

    *weightiwls = *weight * el;

    }

  *tildey =  (*response - el)/el;

  }



void DISTRIBUTION_nbinomial::outoptions(void)
{
    DISTRIBUTION::outoptions();
    ST::string dist;
    if(ver == nb) dist = "negative binomial";
    else if(ver == poga) dist = "poisson-gamma";
    else dist = "poisson-inverse gaussian";
    ST::string proposal;
    if(pscale == gam) proposal = "gamma";
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

void DISTRIBUTION_nbinomial::outresults(void)
{

    DISTRIBUTION::outresults();
    optionsp->out("\n\n");

    double sum=0.0;

    if(ver!=nb)
    {
        if(ver==poig)
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



//Outfile for the nu-sample
//Write out nu

    unsigned i, l = pathresultsscale.length();

    if(ver != nb)
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





    if(ver==poig)
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





void DISTRIBUTION_nbinomial::update(void)
{

    // For hierarchical models, initialize the multiplicative random effects in
    // the first iteration with exp(intercept).

    if(optionsp->get_nriter()==1)
    {
         unsigned i;
         double *nuwork = nu.getV();
         if(hierarchical)
         {
            double *intwork = hierint.getV();
            double expint = exp(*intwork);
            for(i=0; i<nrobs; i++, nuwork++)
            {
                *nuwork = expint;
            }
         }
    }

    // Update of the multiplicative random effects in the POGA and POIG models
    // If hierarchical=True, then update of the hierarchical intercept.


    if(ver!=nb)
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

    // Update of the scale parameter for all the models. Differentiations
    // for the models are included in the function.

    update_scale();

    // Update of the hyperparameter b of the prior for the scale parameter.
    // It is the same function for all the models.

    b_pri(0,0) = update_b_pri();
    b_pri_save.update();

    // Calculate the acceptance rates with the help of the vector acceptance
    // when the running is ended.

    if(optionsp->get_nriter()==optionsp->get_iterations())
    {

        if(ver!=poig)
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


bool DISTRIBUTION_nbinomial::posteriormode(void)
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


bool DISTRIBUTION_nbinomial::posteriormode_converged(const unsigned & itnr)
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



double DISTRIBUTION_nbinomial::update_nu(void)
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




    //Update nu

    // For the Poisson-Gamma Model: Gibbsampling!! Formula (5.15) from diss

    if(ver==poga)
    {

    // nuaux = old nu!

        for (i=0;i<nrobs;i++,nuwork++, respwork++,linwork++, nusavep++)
        {
            mu = exp(*linwork)/(*nuwork);
            nuaux = *nuwork;
            if(hierarchical)
            {
                *nuwork = randnumbers::rand_gamma(*scalework +
                            *respwork, *scalework/expint + mu);
            }
            else *nuwork = randnumbers::rand_gamma(*scalework +
                            *respwork, *scalework + mu);

            *linwork += log((*nuwork)/nuaux);

            *nusavep = *nuwork;
            // For KFZ_Data!!!!
            if(oversize)
            {
                if(i==help && j<10)
                {
                    *nusavekfzp=*nuwork;
                    help += helpsum;
                    nusavekfzp++;
                    j++;
                }
            }
            *sum += *nuwork;
            *sum2 += log(*nuwork);
        }
    }
    // For the Poisson-Inverse Gaussian Model: M-H-Algorithm.
    // Formulae: (5.16), (5.31) and (5.32) from diss.

    else
    {
// nuaux = new nu!
        double priori_ratio;
        double log_likelihood_ratio;
        double proposal_ratio;
        double *pwork = pvar.getV();
        double *acceptwork=accept.getV();
        double alpha;
        acceptwork += 1;
        pwork += 1;
        double lognuaux_nu;
        double diffnu_nuaux;
        double nupost; //maximum of the posterior for nu

        for(i=0; i<nrobs; i++, nuwork++,respwork++, linwork++,
                    acceptwork++, pwork++, nusavep++)
        {



            mu = exp(*linwork)/(*nuwork);

            // (5.31)

            if(hierarchical) nupost = ((*respwork-1.5)+sqrt((*respwork-1.5)*(*respwork-1.5) +
                    2.0*(*scalework)* mu*expint+ *scalework*(*scalework)))/
                    (2.0*mu+*scalework/expint);

            else nupost = ((*respwork-1.5)+sqrt((*respwork-1.5)*(*respwork-1.5) +
                    2.0*(*scalework)* mu+ *scalework*(*scalework)))/
                    (2.0*mu+*scalework);

            // (5.32)

            if(nupost > *pwork)
            {

                nuaux = nupost - *pwork + 2*(*pwork)*(randnumbers::uniform());

            }
            else
            {
                nuaux = (nupost + *pwork)*(randnumbers::uniform());

            }
            proposal_ratio = 0.0;

            lognuaux_nu = log(nuaux/ *nuwork); //help var. for the calculations

            diffnu_nuaux = *nuwork - nuaux;    //help var. for the calculations


            if(hierarchical) priori_ratio = -3/2*lognuaux_nu + *scalework * 0.5 *
                            (diffnu_nuaux/expint + expint/(*nuwork)-expint/nuaux);
            else priori_ratio = -3/2*lognuaux_nu + *scalework * 0.5 *
                            (diffnu_nuaux + 1/(*nuwork)-1/nuaux);
            // (5.16)

            log_likelihood_ratio  =  exp(*linwork)-
                            exp(*linwork +lognuaux_nu) +
                            *respwork * lognuaux_nu;
            alpha = log_likelihood_ratio + priori_ratio + proposal_ratio;
            double h =log(randnumbers::uniform( ));
            if(h<=alpha)
            {
                *linwork += lognuaux_nu;
                *nuwork = nuaux;
                *acceptwork=*acceptwork+1;

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

            *sum += *nuwork;
            *sum2 += 1/(*nuwork);

            // Tuning function for the acceptance rates of the M-H algorithm
            // It will be called only in the burnin phase and each 100 iterations.

            if(optionsp->get_nriter() % 100 == 0 &&
            optionsp->get_nriter() <= optionsp->get_burnin())
            {

            pwork_tunin(i+1);

            }
        }
    }

    return 0;

}


double DISTRIBUTION_nbinomial::update_hierint(void)
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


double DISTRIBUTION_nbinomial::update_scale(void) const
{
    double *scalework = scale.getV();

    // The update for the negative binomial and POGA models is a M-H step!

    if(ver!=poig)
    {
        double scaleaux=*scalework;
        double log_likelihood_ratio;
        double priori_ratio;
        double proposal_ratio;
        double alpha;
        double *acceptwork=accept.getV();

        // Proposal value for scale depending on the proposal distribution
        // And proposal_ratio also depending on the proposal distribution
        // See below.

        proposal_ratio = proposal_scale();

        // After this, *scalework stores the new value and scaleaux the old one

        // For Negative binomial model:

        if(ver==nb)  log_likelihood_ratio = log_nbin(*scalework, scaleaux);

        // For POGA model:

        else
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



double DISTRIBUTION_nbinomial::update_b_pri(void)
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


double DISTRIBUTION_nbinomial::proposal_scale(void) const
{


    double *scalework = scale.getV();
    double scaleaux = *scalework;     //alter Wert!!!!!
    double proposal_ratio;
    double *pwork = pvar.getV();


        // For a uniform proposal: pwork stores the width of the uniform proposal.
        // The central point is given by the old value of scale.
        // Formula (5.38) from diss.

       if(pscale == unif)       // if uniform
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
            while(a_neu < exp(-16*log(10.)))
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

double DISTRIBUTION_nbinomial::pwork_tunin(unsigned i) const
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


// Logarithm of the densitiy Gamma(a, b)

double DISTRIBUTION_nbinomial::log_gamma(double & x,double & a,double & b) const
{
    return a*log(b)-lgamma(a)+(a-1)*log(x)-b*x;
}



double DISTRIBUTION_nbinomial::log_gamma_likelihood(double &s,
                                                    double &s_neu) const
{
    double summe = 0.0;
    double *sum = sum_nu.getV();
    double *sum2 = sum2_nu.getV();
    summe = nrobs*(s_neu*log(s_neu)-s*log(s)+lgamma(s)-lgamma(s_neu)) +
          (s - s_neu)*(*sum - *sum2);
    return summe;
}

double DISTRIBUTION_nbinomial::log_gamma_likelihood_hier(double &s,
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



// Logarithm of the Quotient: G_(a,b)(x_prop)/G_(a, b)(x)

double DISTRIBUTION_nbinomial::log_gamma_quot(double & x, double & x_prop,
                                              const double & a,
                                              const double & b) const
{
    return (a-1)*log(x_prop/x)-b*(x_prop-x);
}


//Loggamma-Funktion

double DISTRIBUTION_nbinomial::lgamma(const double & xx) const
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





double DISTRIBUTION_nbinomial::log_nbin(const double & s_neu,
                                        const double & s) const
{
    unsigned i;
    double *respwork = response.getV();
    double *linwork = (*linpred_current).getV();
    double sum = 0.0;
    double hilfe = exp(-300*log(10.));
    if((s > hilfe) & (s_neu > hilfe))
    {
        for(i=0; i<nrobs; i++, respwork++, linwork++)
            {
                sum += lgamma(s_neu+*respwork)-lgamma(s+*respwork)-
                    (*respwork+s_neu)*log(exp(*linwork)+s_neu)+
                    (*respwork+s)*log(exp(*linwork)+s);
            }
        return sum + nrobs*(lgamma(s)-lgamma(s_neu)+s_neu*log(s_neu)-s*log(s));
    }
    else
    {
        for(i=0; i<nrobs; i++, respwork++, linwork++)
            {
                sum += log_gamma_function_quot(s, s_neu, *respwork)-
                    s_neu*log((exp(*linwork)+s_neu)/s_neu)+
                    s*log((exp(*linwork)+s)/s)+
                    *respwork*log((exp(*linwork)+s)/(exp(*linwork)+s_neu));
            }
        return sum;
    }
}



double DISTRIBUTION_nbinomial::log_gamma_function_quot(const double & g,
                                const double & g_prop, const unsigned & y) const
{
    if(y==0) return 0;
    else
    {
        unsigned i;
        double res=log(g_prop/g);
        for(i=1; i<y; i++)
        {
            res += log((g_prop + i)/(g + i));
        }
     return res;
    }
}


void DISTRIBUTION_nbinomial::add_nu(double m) const
{
    double *nuwork = nu.getV();
    unsigned i;
    for(i=0; i<nrobs; i++, nuwork++)
    {
     *nuwork *= m;
    }
}


} // end: namespace MCMC









