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



#include "tvariance2dim.h"

namespace MCMC
{

FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(MCMCoptions * o,
                   FULLCOND_pspline_surf_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw)
           : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
  {
  spatial = false;
  Laplace = false;

  rowwise = rw;
  Kp = p;
  pathresults = pres;
  nu = v;

  if (!rowwise)
    {
    m = sqrt(static_cast<double>(p->get_nrpar()));
    nrpar = (m-1)*(m-1)*2+2*(m-1);

    SparseMatrix Ksp = Kmrflinear(m,m);
    datamatrix de(m*m-1,1,0.0);
    datamatrix ud = datamatrix(m*m-1,m,0.0);
    for(unsigned i=0;i<m*m-2;i++)
      {
      de(i,0) = Ksp(i,i);
      for(unsigned j=0;j<ud.cols();j++)
        {
        if (i+j+1 < m*m)
          ud(i,j) = Ksp(i,i+j+1);
        }
      }
    de(m*m-2,0) = Ksp(m*m-2,m*m-2);
    K11 = envmatdouble(bandmatdouble(de,ud));
    detalt = K11.getLogDet();
    detneu = 0.0;
    nrrows = unsigned(bs/2);
    betakvec = vector<double>(0);
    rowvec = vector<double>(0);
    colvec = vector<double>(0);
    deltapropvec = vector<double>(0);
    }
  else
    {
    m = sqrt(static_cast<double>(p->get_nrpar()));
    nrpar = m*m;
    u = datamatrix (nrpar,1,0);
    }

  setbeta(nrpar,1,1);

  }


FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(MCMCoptions * o,
                   FULLCOND_nonp_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw)
           : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
  {
  spatial = true;
  Laplace = false;

  rowwise = rw;
  Kp_spat = p;
  pathresults = pres;
  nu = v;

//  dgam = DISTRIBUTION_gamma();

    envmatdouble Kenv = Kp_spat->get_K();

    nrpar = 0;
    for(unsigned i=0;i<Kenv.getDim();i++)
      nrpar += Kenv.getDiag(i);
    nrpar = nrpar/2;

    datamatrix help(Kenv.getDim()-1,Kenv.getDim()-1,0);

    for(unsigned i=0;i<help.rows();i++)
      {
      for(unsigned j=0;j<help.cols();j++)
        {
        help(i,j) = Kenv(i,j);
        }
      }

    unsigned k = 0;
    indexmat = statmatrix<int>(nrpar,3,0);

    for(unsigned i=0;i<help.rows()+1;i++)
      {
      for(unsigned j=0;j<help.cols()+1;j++)
        {
        if(Kenv(i,j) != 0 && j>i)
          {
          indexmat(k,0) = i;
          indexmat(k,1) = j;
          k++;
          }
        }
      }

    K11 = envmatdouble(help);
    detalt = K11.getLogDet();
    detneu = 0.0;
    nrrows = bs;
    betakvec = vector<double>(0);
    rowvec = vector<double>(0);
    colvec = vector<double>(0);
    deltapropvec = vector<double>(0);

  setbeta(nrpar,1,1);

  }


FULLCOND_tvariance2dim::FULLCOND_tvariance2dim(const FULLCOND_tvariance2dim & t)
  : FULLCOND(FULLCOND(t))
  {
//  dgam = t.dgam;
  rowwise = t.rowwise;
  Kp = t.Kp;
  Kp_spat = t.Kp_spat;
  spatial = t.spatial;
  Laplace = t.Laplace;
  indexmat = t.indexmat;
  pathresults = t.pathresults;
  nu = t.nu;
  m = t.m;
  u = t.u;

  K11 = t.K11;
  detalt = t.detalt;
  detneu = t.detneu;
  betakvec = t.betakvec;
  rowvec = t.rowvec;
  colvec = t.colvec;
  deltapropvec = t.deltapropvec;
  nrrows = t.nrrows;
  }


const FULLCOND_tvariance2dim &
FULLCOND_tvariance2dim::operator=(const FULLCOND_tvariance2dim & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
//  dgam = t.dgam;
  rowwise = t.rowwise;
  Kp = t.Kp;
  Kp_spat = t.Kp_spat;
  spatial = t.spatial;
  Laplace = t.Laplace;
  indexmat = t.indexmat;
  pathresults = t.pathresults;
  nu = t.nu;
  m = t.m;
  u = t.u;

  K11 = t.K11;
  detalt = t.detalt;
  detneu = t.detneu;
  betakvec = t.betakvec;
  rowvec = t.rowvec;
  colvec = t.colvec;
  deltapropvec = t.deltapropvec;
  nrrows = t.nrrows;
  return *this;
  }


void FULLCOND_tvariance2dim::update(void)
  {

  unsigned step = 2;

  if(spatial)
    {
    if(optionsp->get_nriter()<optionsp->get_burnin() || (optionsp->get_nriter()-1)%step==0)
      {
      if(Laplace)
        update_spat_laplace();
      else
        update_spat();
      }
    }
  else
    {
    update_2dim();
    }

  }

void FULLCOND_tvariance2dim::update_2dim(void)
  {

  if (!rowwise)
    {

    unsigned i,j,l;
    int k = 0;
    double aneu = double(nu)/2;
    double bneu;

    double alpha,u,betak;
    double deltaprop;
    unsigned row,col;

    unsigned dim = m*m;

    for (row=0;row<dim;row++)
      {

      i = row/m;
      j = row%m;

      col=row+1;
      if(j < m-1)
        {
        betak = beta(k,0);

        bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
        deltaprop = randnumbers::rand_gamma(aneu,bneu);

        deltapropvec.push_back(deltaprop);
        rowvec.push_back(row);
        colvec.push_back(col);
        betakvec.push_back(betak);

        K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
        if(col<K11.getDim())
          {
          K11.set(row,col,-deltaprop);
          K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
          }

        k++;
        }

      col=row+m;
      if(i < m-1)
        {
        betak = beta(k,0);

        bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
        deltaprop = randnumbers::rand_gamma(aneu,bneu);

        deltapropvec.push_back(deltaprop);
        rowvec.push_back(row);
        colvec.push_back(col);
        betakvec.push_back(betak);

        K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
        if(col<K11.getDim())
          {
          K11.set(row,col,-deltaprop);
          K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
          }

        k++;
        }

      if((row+1)%nrrows == 0 || (row+1)==dim)
        {

        if(detalt==detneu)
          K11.decomp2(row-nrrows+1);
        detneu = K11.getLogDet();

        alpha = 0.5*(detneu - detalt);
        u = log(uniform());

        nrtrials++;

        if(u <= alpha)
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
    	    beta(k-deltapropvec.size()+l,0) = deltapropvec[l];
            Kp->setK(rowvec[l],colvec[l],-deltapropvec[l]);
	        }
          detalt = detneu;
          acceptance++;
          }
        else
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
            K11.setDiag(rowvec[l],K11(rowvec[l],rowvec[l]) - deltapropvec[l] + betakvec[l]);
            if(colvec[l]<K11.getDim())
              {
              K11.set(rowvec[l],colvec[l],-betakvec[l]);
              K11.setDiag(colvec[l],K11(colvec[l],colvec[l]) - deltapropvec[l] + betakvec[l]);
              }
    	    }
          }

        deltapropvec = vector<double>(0);
	    rowvec = vector<double>(0);
        colvec = vector<double>(0);
        betakvec = vector<double>(0);

        } // END:       if(row%nrrows == 0)
      } // END:    for (row=0;row<dim;row++)

/*
    unsigned i,j;
    int k = 0;
//    double aneu = double(nu)/2;
    double aneu = double(nu)/2+0.5;
    double bneu;

    double deltaprop,deltaprop2,alpha,u;
    unsigned row,col;

    for (i=0;i<m;i++)
      {
      for(j=0;j<m;j++)
        {

        if ( (i < m-1) && (j < m-1) )
          {

          row = i*m+j;
          col = i*m+j+1;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          deltaprop = randnumbers::rand_gamma(aneu,bneu);

          K11.set(row,col,-deltaprop);
          K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
          K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));

          k++;

          row = i*m+j;
          col = (i+1)*m+j;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          deltaprop2 = randnumbers::rand_gamma(aneu,bneu);

          K11.set(row,col,-deltaprop2);
          K11.setDiag(row,K11(row,row) + deltaprop2 - beta(k,0));
          K11.setDiag(col,K11(col,col) + deltaprop2 - beta(k,0));

          if(detalt==detneu)
            K11.decomp2(row);
          detneu = K11.getLogDet();

          alpha = 0.5*(detneu - detalt);
          alpha += 0.5*( log(beta(k-1,0)) + log(beta(k,0)) - log(deltaprop) - log(deltaprop2) );
          u = log(uniform());
          nrtrials++;

          if(u <= alpha)
            {
            beta(k-1,0) = deltaprop;
            beta(k,0) = deltaprop2;

            Kp->setK(row,col-m+1,-deltaprop);
            Kp->setK(row,col,-deltaprop2);

            detalt = detneu;
            acceptance++;
            }
          else
            {
            col = i*m+j+1;
            K11.set(row,col,-beta(k-1,0));
            K11.setDiag(row,K11(row,row) - deltaprop + beta(k-1,0));
            K11.setDiag(col,K11(col,col) - deltaprop + beta(k-1,0));

            col = (i+1)*m+j;
            K11.set(row,col,-beta(k,0));
            K11.setDiag(row,K11(row,row) - deltaprop2 + beta(k,0));
            K11.setDiag(col,K11(col,col) - deltaprop2 + beta(k,0));
            }

          k++;

          }  // end: if ( (i < m-1) && (j < m-1) )

        if ( (i< m-1) && (j == m-1) )
          {

          row = i*m+j;
          col = (i+1)*m+j;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          deltaprop = randnumbers::rand_gamma(aneu,bneu);

          K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
          if(col<K11.getDim())
            {
            K11.set(row,col,-deltaprop);
            K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
            }

          if(detalt==detneu)
            K11.decomp2(row);
          detneu = K11.getLogDet();

          alpha = 0.5*(detneu - detalt);
          alpha += 0.5*( log(beta(k,0)) - log(deltaprop) );
          u = log(uniform());
          nrtrials++;

          if(u <= alpha)
            {
            beta(k,0) = deltaprop;

            Kp->setK(row,col,-deltaprop);

            detalt = detneu;
            acceptance++;
            }
          else
            {
            K11.setDiag(row,K11(row,row) - deltaprop + beta(k,0));
            if(col<K11.getDim())
              {
              K11.set(row,col,-beta(k,0));
              K11.setDiag(col,K11(col,col) - deltaprop + beta(k,0));
              }
            }

          k++;

          } // end: if ( (i< m-1) && (j == m-1) )

        if ( (i == m-1) && (j < m-1) )
          {

          row = i*m+j;
          col = i*m+j+1;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          deltaprop = randnumbers::rand_gamma(aneu,bneu);

          K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
          if(col<K11.getDim())
            {
            K11.set(row,col,-deltaprop);
            K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
            }

          if(detalt==detneu)
            K11.decomp2(row);
          detneu = K11.getLogDet();

          alpha = 0.5*(detneu - detalt);
          alpha += 0.5*( log(beta(k,0)) - log(deltaprop) );
          u = log(uniform());
          nrtrials++;

          if(u <= alpha)
            {
            beta(k,0) = deltaprop;

            Kp->setK(row,col,-deltaprop);

            detalt = detneu;
            acceptance++;
            }
          else
            {
            K11.setDiag(row,K11(row,row) - deltaprop + beta(k,0));
            if(col<K11.getDim())
              {
              K11.set(row,col,-beta(k,0));
              K11.setDiag(col,K11(col,col) - deltaprop + beta(k,0));
              }
            }

          k++;

          } // end: if ( (i == m-1) && (j < m-1) )

        }  // end: for(j=0;j<m;j++)

      } // end: for (i=0;i<m;i++)
*/
/*
    for (i=0;i<m;i++)
      {
      for(j=0;j<m;j++)
        {

        if ( (i < m-1) && (j < m-1) )
          {

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          *Kmatupper = -beta(k,0);
          k++;

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          Kmatupper+= m-1;
          *Kmatupper = -beta(k,0);
          *Kmatupper++;
          k++;

          }  // end: if ( (i < m-1) && (j < m-1) )

        if ( (i< m-1) && (j == m-1) )
          {

          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i+1,j,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          Kmatupper+= m-1;
          *Kmatupper = -beta(k,0);
          Kmatupper++;
          k++;

          } // end: if ( (i< m-1) && (j == m-1) )

        if ( (i == m-1) && (j < m-1) )
          {
          bneu = nu*0.5 + 0.5*Kp->compute_squareddiff(i,j,i,j+1,m);
          beta(k,0) = randnumbers::rand_gamma(aneu,bneu);
          *Kmatupper = -beta(k,0);
          Kmatupper+= m;
          k++;

          } // end: if ( (i == m-1) && (j < m-1) )

        }  // end: for(j=0;j<m;j++)

      } // end: for (i=0;i<m;i++)
*/

    unsigned size = m*m;
    Kmatdiag = Kp->getdiagpointer();
    Kmatupper = Kp->getupperpointer();
    double * phelp;

    for (i=0;i<size;i++,Kmatdiag++)
      {

      *Kmatdiag = 0;

      if (i < size-1)
        {
        *Kmatdiag += - *Kmatupper;
        Kmatupper += m-1;
        *Kmatdiag += - *Kmatupper;
        Kmatupper ++;
        }

      if (i > 0)
        {

        phelp =  Kp->getupperpointer()+(i-1)*m;

        *Kmatdiag += - *phelp;

        if (i >= m)
          {

          phelp =  Kp->getupperpointer()+(i-m)*m + m-1;
          *Kmatdiag += - *phelp;

          }

        } // end: if (i > 0)

      }

    }         // end: !rowwise
  else        // rowwise
    {

    unsigned i,j;

    Kp->compute_squareddiff(u);

    Kmatdiag = Kp->getdiagpointer();
    Kmatupper = Kp->getupperpointer();
    double * workbeta = beta.getV()+1;
    double * worku = u.getV()+1;
    double wold,wnew;
    double v = nu/2.0;

    for (i=0;i<m;i++,workbeta++,worku++,Kmatdiag++,Kmatupper++)
      {
      *workbeta = rand_invgamma(v+0.5,v+0.5* *worku);
      wold=1.0/(*workbeta);
      *Kmatdiag = wold;
      *Kmatupper = -wold;
      Kmatdiag++;
      Kmatupper++;
      workbeta++;
      worku++;

      for (j=1;j<m-1;j++,Kmatdiag++,Kmatupper++,workbeta++,worku++)
        {
        *workbeta = rand_invgamma(v+0.5,v+0.5* *worku);
        wnew = 1.0/(*workbeta);
        *Kmatdiag = wold+wnew;
        *Kmatupper = -wnew;
        wold = wnew;
        }

      *Kmatdiag = wold;
      *Kmatupper = 0;


      }

    acceptance++;

    }

  FULLCOND::update();
  }


void FULLCOND_tvariance2dim::update_spat(void)
  {

    unsigned l;
    int k = 0;
    double aneu = 0.5*double(nu);
    double bneu;

    double alpha,u,betak;
    double deltaprop;
    unsigned row,col;

    while(k<nrpar)
      {

      betak = beta(k,0);
      row = indexmat(k,0);
      col = indexmat(k,1);

      bneu = nu*0.5 + 0.5*Kp_spat->compute_squareddiff(row,col);
      deltaprop = randnumbers::rand_gamma(aneu,bneu);

      deltapropvec.push_back(deltaprop);
      rowvec.push_back(row);
      colvec.push_back(col);
      betakvec.push_back(betak);

      K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
      if(col<K11.getDim())
        {
        K11.set(row,col,-deltaprop);
        K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
        }

      k++;

      if((row+1)%nrrows == 0 || k==nrpar)
        {

        if(detalt==detneu)
          K11.decomp2(row-nrrows+1);
        detneu = K11.getLogDet();

        alpha = 0.5*(detneu - detalt);
        u = log(uniform());

        nrtrials++;

        if(u <= alpha)
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
    	    beta(k-deltapropvec.size()+l,0) = deltapropvec[l];
            Kp_spat->setK(rowvec[l],colvec[l],-deltapropvec[l]);

            Kp_spat->setK(colvec[l],colvec[l],Kp_spat->get(colvec[l],colvec[l])+deltapropvec[l]-betakvec[l]);
            Kp_spat->setK(rowvec[l],rowvec[l],Kp_spat->get(rowvec[l],rowvec[l])+deltapropvec[l]-betakvec[l]);

	        }
          detalt = detneu;
          acceptance++;
          }
        else
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
            K11.setDiag(rowvec[l],K11(rowvec[l],rowvec[l]) - deltapropvec[l] + betakvec[l]);
            if(colvec[l]<K11.getDim())
              {
              K11.set(rowvec[l],colvec[l],-betakvec[l]);
              K11.setDiag(colvec[l],K11(colvec[l],colvec[l]) - deltapropvec[l] + betakvec[l]);
              }
    	    }
          }

        deltapropvec = vector<double>(0);
	    rowvec = vector<double>(0);
        colvec = vector<double>(0);
        betakvec = vector<double>(0);

        } // END:       if(row%nrrows == 0)
      } // END:    for (row=0;row<dim;row++)
/*
ofstream out("c:\\bayesx\\K11.raw");
K11.print4(out);
out.close();

ofstream out2("c:\\bayesx\\Kenv.raw");
(Kp_spat->get_K()).print4(out2);
out2.close();
*/

/*/ Diagonalelemente ausrechnen

  unsigned size = (Kp_spat->get_K()).getDim();

  double help;
  k=0;
  for(unsigned i=0;i<size;i++)
    {
    help = 0.0;
    row = indexmat(i,0);
    while(indexmat(i,0) == row)
      {
      help += beta(k,0);
      i++;
      }
    Kp_spat->setK(row,row,-help);
    k++;
    }
*/
//  acceptance++;

  FULLCOND::update();
  }


void FULLCOND_tvariance2dim::update_spat_laplace(void)
  {
/*         geht nicht !!!
    unsigned l;
    int k = 0;
    double aneu = 1.0 + 0.5*double(nu);
    double bneu;

    double alpha,u,betak;
    double deltaprop;
    unsigned row,col;

    double quadform;
    double nu_K = fabs(0.5*(2-double(nrpar)));
    double propnew=0.0;
    double propold=0.0;
    double lognew=0.0;
    double logold=0.0;

//    double proposalvar = 0.15;

    while(k<nrpar)
      {

      betak = beta(k,0);
      row = indexmat(k,0);
      col = indexmat(k,1);

      bneu = 0.5*nu + Kp_spat->compute_fabsdiff(row,col);
      deltaprop = randnumbers::rand_gamma(aneu,bneu);
      lognew  += (0.5*nu-1)*log(deltaprop) - 0.5*nu*deltaprop;
      logold  += (0.5*nu-1)*log(betak)     - 0.5*nu*betak;

      propold +=   (aneu-1)*log(deltaprop) -   bneu*deltaprop;
//      propold += aneu*log(bneu) - dgam.lgammafunc(aneu);

      propnew +=   (aneu-1)*log(betak)     -   bneu*betak;
//      propnew += aneu*log(bneu) - dgam.lgammafunc(aneu);

      deltapropvec.push_back(deltaprop);
      rowvec.push_back(row);
      colvec.push_back(col);
      betakvec.push_back(betak);

      K11.setDiag(row,K11(row,row) + deltaprop - beta(k,0));
      if(col<K11.getDim())
        {
        K11.set(row,col,-deltaprop);
        K11.setDiag(col,K11(col,col) + deltaprop - beta(k,0));
        }

      k++;

      if((row+1)%nrrows == 0 || k==nrpar)
        {

        if(detalt==detneu)
          K11.decomp2(row-nrrows+1);
        detneu = K11.getLogDet();

        quadform = Kp_spat->compute_quadform();
        if(quadform>0.1)
          logold += 0.5*nu_K*log(quadform);
        if(nu_K<100)
          logold += log(besselK(sqrt(2*quadform),nu_K));
        else
          logold += log_besselK(sqrt(2*quadform),nu_K);

        for(l=0;l<deltapropvec.size();l++)
          {
          Kp_spat->setK(rowvec[l],colvec[l],-deltapropvec[l]);

          Kp_spat->setK(colvec[l],colvec[l],Kp_spat->get(colvec[l],colvec[l])+deltapropvec[l]-betakvec[l]);
          Kp_spat->setK(rowvec[l],rowvec[l],Kp_spat->get(rowvec[l],rowvec[l])+deltapropvec[l]-betakvec[l]);
          }

        quadform = Kp_spat->compute_quadform();
        if(quadform>0.1)
          lognew += 0.5*nu_K*log(quadform);
        if(nu_K<100)
          lognew += log(besselK(sqrt(2*quadform),nu_K));
        else
          lognew += log_besselK(sqrt(2*quadform),nu_K);

        alpha = 0.5*(detneu - detalt) + lognew + propnew - logold - propold;
        u = log(uniform());

        nrtrials++;

        if(u <= alpha)
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
    	    beta(k-deltapropvec.size()+l,0) = deltapropvec[l];
//            Kp_spat->setK(rowvec[l],colvec[l],-deltapropvec[l]);

//            Kp_spat->setK(colvec[l],colvec[l],Kp_spat->get(colvec[l],colvec[l])+deltapropvec[l]-betakvec[l]);
//            Kp_spat->setK(rowvec[l],rowvec[l],Kp_spat->get(rowvec[l],rowvec[l])+deltapropvec[l]-betakvec[l]);

	        }
          detalt = detneu;
          acceptance++;
          }
        else
          {
          for(l=0;l<deltapropvec.size();l++)
	        {
            Kp_spat->setK(rowvec[l],colvec[l],-betakvec[l]);
            Kp_spat->setK(colvec[l],colvec[l],Kp_spat->get(colvec[l],colvec[l]) - deltapropvec[l] + betakvec[l]);
            Kp_spat->setK(rowvec[l],rowvec[l],Kp_spat->get(rowvec[l],rowvec[l]) - deltapropvec[l] + betakvec[l]);

            K11.setDiag(rowvec[l],K11(rowvec[l],rowvec[l]) - deltapropvec[l] + betakvec[l]);
            if(colvec[l]<K11.getDim())
              {
              K11.set(rowvec[l],colvec[l],-betakvec[l]);
              K11.setDiag(colvec[l],K11(colvec[l],colvec[l]) - deltapropvec[l] + betakvec[l]);
              }
    	    }
          }

        lognew = 0.0;
        logold = 0.0;
        propnew = 0.0;
        propold = 0.0;

        deltapropvec = vector<double>(0);
	    rowvec = vector<double>(0);
        colvec = vector<double>(0);
        betakvec = vector<double>(0);

        } // END:       if(row%nrrows == 0)
      } // END:    for (row=0;row<dim;row++)
*/
  Kp_spat->set_delta(beta);

  FULLCOND::update();
  }


void FULLCOND_tvariance2dim::outresults(void)
  {

  FULLCOND::outresults();

  optionsp->out("  Results are stored in file " + pathresults + "\n");
  optionsp->out("\n");

  unsigned i;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  ofstream outres(pathresults.strtochar());

  ST::string name = title;

  outres << "intnr" << "   ";
  outres << name << "   ";
  outres << "pmean   ";
  outres << "pqu"  << l1  << "   ";
  outres << "pqu"  << l2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << u1  << "   ";
  outres << "pqu"  << u2  << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  double * workbetaqu50 = betaqu50.getV();

  for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
                           workbetaqu_l2_lower_p++,workbetaqu50++,
                           workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
    {
    outres << (i+1) << "   ";
    outres << (i+1) << "   ";
    outres << *workmean << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workbetaqu50 << "   ";
    outres << *workbetaqu_l2_upper_p << "   ";
    outres << *workbetaqu_l1_upper_p << "   ";
    outres << endl;
    }

  }


void FULLCOND_tvariance2dim::outoptions(void)
  {
  optionsp->out("  OPTIONS FOR NONPARAMETRIC TERM: " + title +
                " (weights)\n",true);
  optionsp->out("\n");

  optionsp->out("  Hyperprior nu for variance parameter: " +
                ST::inttostring(nu) + "\n" );
//  optionsp->out("  Blocksize for updating variances: ");
  if(spatial)
    optionsp->out("  Blocksize for updating variances: " + ST::inttostring(nrrows) + " rows of penalty matrix\n" );
  else
    optionsp->out("  Blocksize for updating variances: " + ST::inttostring(nrrows*2) + "\n" );
  optionsp->out("\n");

  }


double besselK(const double x, const double xnu)
  {
  double x1,x2,x3,x4;
  bessik(x,xnu,x1,x2,x3,x4);
  return x2;
  }

void bessik(const double x, const double xnu, double &ri, double &rk, double &rip, double &rkp)
  {
  const int MAXIT=10000;
  const double EPS=std::numeric_limits<double>::epsilon();
  const double FPMIN=std::numeric_limits<double>::min()/EPS;
  const double XMIN=2.0;
//  const double PI=3.141592653589793;
  double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
          gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
          ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;
  int i,l,nl;

  assert(x > 0.0 && xnu >= 0.0);
  nl=int(xnu+0.5);
  xmu=xnu-nl;
  xmu2=xmu*xmu;
  xi=1.0/x;
  xi2=2.0*xi;
  h=xnu*xi;
  if (h < FPMIN) h=FPMIN;
  b=xi2*xnu;
  d=0.0;
  c=h;
  for (i=0;i<MAXIT;i++) {
    b += xi2;
    d=1.0/(b+d);
    c=b+1.0/c;
    del=c*d;
    h=del*h;
    if (fabs(del-1.0) <= EPS) break;
  }
  assert(i<MAXIT);
  ril=FPMIN;
  ripl=h*ril;
  ril1=ril;
  rip1=ripl;
  fact=xnu*xi;
  for (l=nl-1;l >= 0;l--) {
    ritemp=fact*ril+ripl;
    fact -= xi;
    ripl=fact*ritemp+ril;
    ril=ritemp;
  }
  f=ripl/ril;
  if (x < XMIN) {
    x2=0.5*x;
    pimu=PI*xmu;
    fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
    d = -log(x2);
    e=xmu*d;
    fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
    beschb(xmu,gam1,gam2,gampl,gammi);
    ff=fact*(gam1*cosh(e)+gam2*fact2*d);
    sum=ff;
    e=exp(e);
    p=0.5*e/gampl;
    q=0.5/(e*gammi);
    c=1.0;
    d=x2*x2;
    sum1=p;
    for (i=1;i<=MAXIT;i++) {
      ff=(i*ff+p+q)/(i*i-xmu2);
      c *= (d/i);
      p /= (i-xmu);
      q /= (i+xmu);
      del=c*ff;
      sum += del;
      del1=c*(p-i*ff);
      sum1 += del1;
      if (fabs(del) < fabs(sum)*EPS) break;
    }
    assert(i <= MAXIT);
    rkmu=sum;
    rk1=sum1*xi2;
  } else {
    b=2.0*(1.0+x);
    d=1.0/b;
    h=delh=d;
    q1=0.0;
    q2=1.0;
    a1=0.25-xmu2;
    q=c=a1;
    a = -a1;
    s=1.0+q*delh;
    for (i=1;i<MAXIT;i++) {
      a -= 2*i;
      c = -a*c/(i+1.0);
      qnew=(q1-b*q2)/a;
      q1=q2;
      q2=qnew;
      q += c*qnew;
      b += 2.0;
      d=1.0/(b+a*d);
      delh=(b*d-1.0)*delh;
      h += delh;
      dels=q*delh;
      s += dels;
      if (fabs(dels/s) <= EPS) break;
    }
    assert(i < MAXIT);
    h=a1*h;
    rkmu=sqrt(PI/(2.0*x))*exp(-x)/s;
    rk1=rkmu*(xmu+x+0.5-h)*xi;
  }
  rkmup=xmu*xi*rkmu+rk1;
  rimu=xi/(f*rkmu-rkmup);
  ri=(rimu*ril1)/ril;
  rip=(rimu*rip1)/ril;
  for (i=1;i <= nl;i++) {
    rktemp=(xmu+i)*xi2*rk1+rkmu;
    rkmu=rk1;
    rk1=rktemp;
  }
  rk=rkmu;
  rkp=xnu*xi*rkmu+rk1;
  }

void beschb(const double x, double &gam1, double &gam2, double &gampl, double &gammi)
  {

  const int NUSE1=7, NUSE2=8;
  static const double c1_d[7] = {
      -1.142022680371168e0,6.5165112670737e-3,
      3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
      3.67795e-11,-1.356e-13};
  static const double c2_d[8] = {
      1.843740587300905e0,-7.68528408447867e-2,
      1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
      2.423096e-10,-1.702e-13,-1.49e-15};
  double xx;
  static Vec_DP c1(c1_d,7),c2(c2_d,8);

  xx=8.0*x*x-1.0;
  gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
  gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
  gampl= gam2-x*gam1;
  gammi= gam2+x*gam1;

  }

double chebev(const double a, const double b, Vec_I_DP &c, const int m, const double x)
  {
  double d=0.0,dd=0.0,sv,y,y2;
  int j;

  assert((x-a)*(x-b) <= 0.0);
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for (j=m-1;j>0;j--) {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];

  }

double log_besselK(const double x, const double xnu)
  {
// nach Abramowitz/Stegun 9.7.8
  double z,eta,u1,u2,t,t2;

  z   = x/xnu;
  t   = 1.0/sqrt(1.0+z*z);
  t2  = t*t;
  eta = sqrt(1.0+z*z)+log(z/(1+sqrt(1.0+z*z)));
  u1  = t*(3.0-5.0*t2)/24.0;
  u2  = t2*(81.0-462.0*t2+385.0*t2*t2)/1152.0;

  z = -xnu*eta - 0.25*log(1+z*z) + log(1.0 - u1/xnu + u2/(xnu*xnu));

  return  z;
  }


} // end: namespace MCMC




