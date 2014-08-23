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



#include "design_kriging.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//------------ CLASS: DESIGN_kriging implementation of member functions --------
//------------------------------------------------------------------------------


void DESIGN_kriging::compute_orthogonaldecomp(void)
  {

  compute_XtransposedWX();
  datamatrix R = XWXfull.root();

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\XWX.res");
  // R.prettyPrint(out);
  // TEST

  datamatrix Rt = R.transposed();
  datamatrix RinvKRtinv = R.inverse()*Kfull*Rt.inverse();

  s = datamatrix(nrpar,1,0);

  bool ecorrect = eigen2(RinvKRtinv,s);
  eigensort(s,RinvKRtinv);

  QtRinv = RinvKRtinv.transposed()*R.inverse();
  RtinvQ = Rt.inverse()*RinvKRtinv;

  u = datamatrix(nrpar,1,0);

  }


void DESIGN_kriging::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  */

  datanames = vn;

  int f;

  f = op[2].strtolong(nrknots);
  if (op[22]=="0.5")
    nu = 0.5;
  else if (op[22]=="1.5")
    nu = 1.5;
  else if (op[22]=="2.5")
    nu = 2.5;
  else if (op[22]=="3.5")
    nu = 3.5;

  f = op[23].strtodouble(maxdist);

  if(maxdist<=0) // wähle maxdist so, dass Korrelation für
                 // Punkte mitmaximalem Abstand = 0.0001
    {
    if(nu==0.5)
      {
      maxdist=9.21034037;//4.605170186;
      }
    else if(nu==1.5)
      {
      maxdist=11.75637122;//6.638352068;
      }
    else if(nu==2.5)
      {
      maxdist=13.53592464;//8.022007057;
      }
    else if(nu==3.5)
      {
      maxdist=15.01510426;//9.158140446;
      }
    }

  knotdatapath = op[36];

  }


DESIGN_kriging::DESIGN_kriging(void) : DESIGN()
  {

  }


void DESIGN_kriging::read_knots_from_data(void)
  {
  if (knotdatapath!="")
    {
    datamatrix help;
    ifstream in(knotdatapath.strtochar());
    if(!in.fail())
      {
      help.prettyScan(in);
      nrknots=help.rows();
      unsigned i;
      for(i=0; i<nrknots; i++)
        {
        xknots.push_back(help(i,0));
        yknots.push_back(help(i,1));
        }
      }

    }
  }


void DESIGN_kriging::help_construct(const datamatrix & dmy, const datamatrix & iv,
                               vector<ST::string> & op, vector<ST::string> & vn,
                               datamatrix & kd)
  {

  center = false;

  full = true;

  read_options(op,vn);

  type = Grf;

  init_data(dmy,iv);

  if (kd.rows() > 1)
    {
    unsigned i;
    nrknots = kd.rows();
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(kd(i,0));
      yknots.push_back(kd(i,1));
      }
    }
  else
    compute_knots(xvalues,yvalues,nrknots,-20,20,xknots,yknots);

  nrpar = xknots.size();

  if (knotdatapath.isvalidfile() != 1)
    {
    ofstream out(knotdatapath.strtochar());
    unsigned i;
    out << datanames[0] << "   " << datanames[1] << endl;
    for(i=0; i<nrknots; i++)
      {
      out << xknots[i] << "   " << yknots[i] << endl;
      }
    out.close();
    }

  // berechne rho
  unsigned i,j;
  rho=0;
  double norm2;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      norm2=(xvalues[i]-xvalues[j])*(xvalues[i]-xvalues[j])+
      (yvalues[i]-yvalues[j])*(yvalues[i]-yvalues[j]);
      if(norm2>rho)
        {
        rho=norm2;
        }
      }
    }
  rho=sqrt(rho)/maxdist;

  compute_penalty();

  compute_tildeZ();

  XWres = datamatrix(nrpar,1);

  XWXfull = datamatrix(nrpar,nrpar);
  WsumtildeZ = datamatrix(Zout.rows(),Zout.cols());
  Wsum = datamatrix(posbeg.size(),1,1);


  }


DESIGN_kriging::DESIGN_kriging(const datamatrix & dm,const datamatrix & iv,
                               const MAP::map & mp,
                               GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                               vector<ST::string> & op, vector<ST::string> & vn)
                               : DESIGN(o,dp,fcl)
  {

  ma = mp;
  mapexisting = true;
  ma.sortmap();
  regions=datamatrix(ma.get_nrregions(),3,0);
  unsigned i;
  ST::string rn;
  MAP::region rh;
  for (i=0;i<regions.rows();i++)
    {
    rn = ma.getname(i);
    rh = ma.get_region(i);
    rn.strtodouble(regions(i,0));
    regions(i,1) = rh.get_xcenter();
    regions(i,2) = rh.get_ycenter();
    }

  // test
  // ofstream out("c:\\bayesx\\testh\\results\\regions.res");
  // regions.prettyPrint(out);
  // ende test


  statmatrix<int> ind(dm.rows(),1);
  ind.indexinit();
  dm.indexsort(ind,0,dm.rows()-1,0,0);

  datamatrix dmy(dm.rows(),2,0);
  vector<ST::string> evh;
  unsigned j=0;
  evh.push_back(ST::doubletostring(regions(j,0)));
  bool found;
  for (i=0;i<dm.rows();i++)
    {
    found=false;
    while (found==false)
      {
      if (dm(ind(i,0),0) == regions(j,0))
        {
        dmy(ind(i,0),0) = regions(j,1);
        dmy(ind(i,0),1) = regions(j,2);
        found=true;
        }
      else
        {
        j++;
        evh.push_back(ST::doubletostring(regions(j,0)));
        }
      }
    }

  // test
  // ofstream out("c:\\bayesx\\testh\\results\\dmy.res");
  // dmy.prettyPrint(out);
  // ende test


  datamatrix kd;
  help_construct(dmy,iv,op,vn,kd);



  ST::string h;
  for(i=0;i<effectvalues.size();i++)
   {
   h =   ST::doubletostring(dm(index_data(posbeg[i],0),0));
   effectvalues[i] = h + "   " + effectvalues[i];
   }

  datanames[0] = datanames[0] + "   xcoord    ycoord";

  }


DESIGN_kriging::DESIGN_kriging(const datamatrix & dm,const datamatrix & iv,
                               GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                               vector<ST::string> & op, vector<ST::string> & vn,
                               datamatrix & kd)
                               : DESIGN(o,dp,fcl)
  {

  help_construct(dm,iv,op,vn,kd);

  /*
  center = false;
  full = true;

  read_options(op,vn);

  type = Grf;

  init_data(dm,iv);

  compute_knots(xvalues,yvalues,nrknots,-20,20,xknots,yknots);
  //  xknots = xvalues;
  //  yknots = yvalues;
  nrpar = xknots.size();

  // berechne rho
  unsigned i,j;
  rho=0;
  double norm2;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      norm2=(xvalues[i]-xvalues[j])*(xvalues[i]-xvalues[j])+
      (yvalues[i]-yvalues[j])*(yvalues[i]-yvalues[j]);
      if(norm2>rho)
        {
        rho=norm2;
        }
      }
    }
  rho=sqrt(rho)/maxdist;

  compute_penalty();

  compute_tildeZ();

  XWres = datamatrix(nrpar,1);

  XWXfull = datamatrix(nrpar,nrpar);
  WsumtildeZ = datamatrix(Zout.rows(),Zout.cols());
  Wsum = datamatrix(posbeg.size(),1,1);
  */

  }


  // COPY CONSTRUCTOR

DESIGN_kriging::DESIGN_kriging(const DESIGN_kriging & m)
    : DESIGN(DESIGN(m))
  {

  ma = m.ma;
  mapexisting = m.mapexisting;
  regions = m.regions;

  knotdatapath = m.knotdatapath;

  xknots = m.xknots;
  yknots = m.yknots;
  xvalues = m.xvalues;
  yvalues = m.yvalues;
  nrknots = m.nrknots;
  rho = m.rho;
  nu = m.nu;
  maxdist = m.maxdist;
  tildeZ_t = m.tildeZ_t;
  XWXfull = m.XWXfull;
  WsumtildeZ = m.WsumtildeZ;
  Kfull = m.Kfull;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_kriging & DESIGN_kriging::operator=(const DESIGN_kriging & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));

  ma = m.ma;
  mapexisting = m.mapexisting;
  regions = m.regions;

  knotdatapath = m.knotdatapath;

  xknots = m.xknots;
  yknots = m.yknots;
  xvalues = m.xvalues;
  yvalues = m.yvalues;
  nrknots = m.nrknots;
  rho = m.rho;
  nu = m.nu;
  maxdist = m.maxdist;
  tildeZ_t = m.tildeZ_t;
  XWXfull = m.XWXfull;
  WsumtildeZ = m.WsumtildeZ;
  Kfull = m.Kfull;
  return *this;
  }


void DESIGN_kriging::compute_knots(const vector<double> & xvals,
                                   const vector<double> & yvals,
                                   unsigned nrknots,double p,double q,
                                   vector<double> & xknots,
                                   vector<double> & yknots)
  {

  int maxsteps=100;
  unsigned nrdiffobs = posbeg.size();

  if(nrknots>xvals.size())
    {
//    errors.push_back("ERROR: More knots requested than different locations observed");
    }
  else if(nrknots==xvals.size())
    {
    xknots=xvals;
    yknots=yvals;
    }
  else
    {
//    optionsp->out("\n");
//    optionsp->out("\n");
//    optionsp->out("Computing knots (this may take some time)\n",true);
//    optionsp->out("\n");
    //Space-Filling-Algorithmus
    unsigned i,j,k;

    unsigned nrcand = nrdiffobs-nrknots; // Anzahl der Kandidatenpunkte

    vector<unsigned> knotindex(nrknots,0);
    vector<unsigned> candindex(nrcand,0);
    datamatrix distmat(nrcand, nrknots,0);     // Distanzmatrix
    datamatrix rowsums(nrcand,1,0);            // Zeilensummen der Distanzmatrix

    datamatrix rowsumswithouti(nrcand,1,0);    // Zeilensummen ohne den Beitrag von Knoten i
    double partialnewrow;   // Zeilensumme der neuen Zeile der Distanzmatrix, die sich ergibt, wenn Knoten i
                            // zum Kandidaten wird (ohne Beitrag des neuen Knotens)

    int swapindex;                        // Index des zu tauschenden Kandidaten
                                          // (Kandidat mit dem minimalen Kriterium beim Loop über den aktuellen Knoten)

    double covcritold=0; // Coverage-Kriterium (ändert sich nur nach Durchlauf über Knoten)
    double covcritnew=-1;   // Neues Coverage-Kriterium (ändert sich nur nach Durchlauf über Knoten)
    double covcritoldi;  // Coverage-Kriterium vorm Swappen von Knoten i (ändert sich nach jedem Swap)
    double covcritnewi;  // Neues Coverage-Kriterium beim Swappen von Knoten i (ändert sich nach jedem Swap)

    // Startdesign und Kandidaten zufällig bestimmen

    datamatrix u(nrdiffobs,1,0);
    statmatrix<int> ind(nrdiffobs,1,0);
    ind.indexinit();
    for(i=0; i<nrdiffobs; i++)
      {
      u(i,0) = uniform();
      }
    u.indexsort(ind,0,nrdiffobs-1,0,0);
    for(i=0; i<nrknots; i++)
      {
      knotindex[i] = ind(i,0);
      }
    for(i=nrknots, j=0; i<nrdiffobs; i++, j++)
      {
      candindex[j] = ind(i,0);
      }

    // Zentriere xvals und yvals
    double xmean=0;
    double ymean=0;
    for(i=0; i<nrdiffobs; i++)
      {
      xmean += xvals[i];
      ymean += yvals[i];
      }
    vector<double> xvalscentered(nrdiffobs);
    vector<double> yvalscentered(nrdiffobs);
    for(i=0; i<nrdiffobs; i++)
      {
      xvalscentered[i] = xvals[i]-xmean;
      yvalscentered[i] = yvals[i]-ymean;
      }

    // Distanzmatrix mit Startdesign und Zeilensummen ausrechnen

    for(i=0; i<nrcand; i++)
      {
      for(j=0; j<nrknots; j++)
        {
        distmat(i,j) = pow(
               (xvalscentered[knotindex[j]]-xvalscentered[candindex[i]])*(xvalscentered[knotindex[j]]-xvalscentered[candindex[i]])
               + (yvalscentered[knotindex[j]]-yvalscentered[candindex[i]])*(yvalscentered[knotindex[j]]-yvalscentered[candindex[i]]),p/2);
        rowsums(i,0) += distmat(i,j);
        }
      }
    for(i=0; i<nrcand; i++)
      {
      covcritold += pow(rowsums(i,0),q/p);
      }
    covcritold = pow(covcritold,1/q);
    covcritoldi = covcritold;

    // Swappen
    unsigned steps=1;
    while(covcritnew<covcritold && steps<=maxsteps)
      {
      if(steps>1)
        {
        covcritold=covcritnew;
        }
      for(i=0; i<nrknots; i++)
        {
        swapindex = -1;
        for(k=0; k<nrcand; k++)
          {
          rowsumswithouti(k,0)=rowsums(k,0)-distmat(k,i);
          }
        partialnewrow=0;
        for(k=0; k<i; k++)
          {
          partialnewrow += pow(
               (xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])*(xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])
               + (yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]])*(yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]]),p/2);
          }
        for(k=i+1; k<nrknots; k++)
          {
          partialnewrow += pow(
               (xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])*(xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])
               + (yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]])*(yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]]),p/2);
          }

        for(j=0; j<nrcand; j++)
          {
          // Loop über die Kandidaten. Jeweil das Coverage-Kriterium bei Swap i vs. j ausrechnen
          // und gegebenenfalls swapindex umsetzen.
          covcritnewi = 0;
          // Loopen über Zeilensummen ohne j
          for(k=0; k<j; k++)
            {
            covcritnewi += pow(rowsumswithouti(k,0) +
               pow((xvalscentered[candindex[k]]-xvalscentered[candindex[j]])*(xvalscentered[candindex[k]]-xvalscentered[candindex[j]])
               + (yvalscentered[candindex[k]]-yvalscentered[candindex[j]])*(yvalscentered[candindex[k]]-yvalscentered[candindex[j]]),p/2),q/p);
            }
          for(k=j+1; k<nrcand; k++)
            {
            covcritnewi += pow(rowsumswithouti(k,0) +
               pow((xvalscentered[candindex[k]]-xvalscentered[candindex[j]])*(xvalscentered[candindex[k]]-xvalscentered[candindex[j]])
               + (yvalscentered[candindex[k]]-yvalscentered[candindex[j]])*(yvalscentered[candindex[k]]-yvalscentered[candindex[j]]),p/2),q/p);
            }
          // Neue Zeile
          covcritnewi += pow(partialnewrow + distmat(j,i),q/p);
          covcritnewi = pow(covcritnewi,1/q);

          if(covcritnewi<covcritoldi)
            {
            swapindex=j;
            covcritnew=covcritnewi;
            covcritoldi=covcritnew;
            }
          }
        if(swapindex!=-1)
          {
          // Swappen
          unsigned h = knotindex[i];
          knotindex[i] = candindex[swapindex];
          candindex[swapindex] = h;

//          unsigned test1=h;
//          unsigned test2=knotindex[i];

          // Distanzmatrix und Zeilensummen updaten
          // Neue Spalte ausrechnen (Index i)
          for(j=0; j<nrcand; j++)
            {
            distmat(j,i) = pow(
               (xvalscentered[knotindex[i]]-xvalscentered[candindex[j]])*(xvalscentered[knotindex[i]]-xvalscentered[candindex[j]])
               + (yvalscentered[knotindex[i]]-yvalscentered[candindex[j]])*(yvalscentered[knotindex[i]]-yvalscentered[candindex[j]]),p/2);
            }

          // Neue Zeile ausrechnen (Index swapindex)
          for(j=0; j<nrknots; j++)
            {
            distmat(swapindex,j) = pow(
               (xvalscentered[knotindex[j]]-xvalscentered[candindex[swapindex]])*(xvalscentered[knotindex[j]]-xvalscentered[candindex[swapindex]])
               + (yvalscentered[knotindex[j]]-yvalscentered[candindex[swapindex]])*(yvalscentered[knotindex[j]]-yvalscentered[candindex[swapindex]]),p/2);
            }

        // Zeilensummen
          for(k=0; k<nrcand; k++)
            {
            rowsums(k,0)=0;
            for(j=0; j<nrknots; j++)
              {
              rowsums(k,0) += distmat(k,j);
              }
            }
          }
        }
      steps++;
      }


    // Knoten speichern
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(xvals[knotindex[i]]);
      yknots.push_back(yvals[knotindex[i]]);
      }
    }
  }




void DESIGN_kriging::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  unsigned j;

  // 1. Indexsort of data
  if (index_data.rows() <= 1)
    {
    index_data = statmatrix<int>(dm.rows(),1);
    index_data.indexinit();
    dm.indexsort2d(index_data,0,dm.rows()-1,0,1,0);
    }

  // TEST
  /*
   ofstream out("c:\\bayesx\\testh\\results\\sortdata.res");
   unsigned t;
   for(t=0;t<dm.rows();t++)
    out << dm(index_data(t,0),0) << "  " << dm(index_data(t,0),1) << endl;
  */
  // TEST


  //2. data = sorted observations, init intvar and intvar2
  data = datamatrix(dm.rows(),2);
  double * workdata = data.getV();
  int * workindex = index_data.getV();
  for (j=0;j<dm.rows();j++,workdata++,workindex++)
    {
    *workdata = dm(*workindex,0);
    workdata++;
    *workdata = dm(*workindex,1);
    }

  if (iv.rows() == dm.rows())
    {
    intvar = iv;
    intvar2 = datamatrix(iv.rows(),1);
    double * workintvar2 = intvar2.getV();
    double * workintvar = intvar.getV();

    for (j=0;j<iv.rows();j++,workintvar++,workintvar2++)
      {
      *workintvar2 = pow(*workintvar,2);
      }
    }

  // 3. Creates posbeg, posend
  posbeg.erase(posbeg.begin(),posbeg.end());
  posend.erase(posend.begin(),posend.end());
  xvalues.erase(xvalues.begin(),xvalues.end());
  yvalues.erase(yvalues.begin(),yvalues.end());
  posbeg.push_back(0);
  double help1 = data(0,0);
  double help2 = data(0,1);
  xvalues.push_back(help1);
  yvalues.push_back(help2);
  for(j=1;j<data.rows();j++)
    {
    if (  data(j,0) != help1 || data(j,1) != help2)
      {
      posend.push_back(j-1);
      if (j < data.rows())
        posbeg.push_back(j);
      xvalues.push_back(data(j,0));
      yvalues.push_back(data(j,1));
      }

    help1 = data(j,0);
    help2 = data(j,1);

    }

  if (posend.size() < posbeg.size())
    posend.push_back(data.rows()-1);

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\xvalues.res");
  // for (j=0;j<xvalues.size();j++)
  //   out << xvalues[j] << "  " << yvalues[j] << endl;

  // ofstream out("c:\\bayesx\\testh\\results\\posbeg.res");
  // for (j=0;j<posbeg.size();j++)
  //   out << posbeg[j] << "  " << posend[j] << endl;
  // TEST

  // 4. initializes ind
  int k;
  workindex = index_data.getV();
  ind = statmatrix<unsigned>(dm.rows(),1);
  for (j=0;j<posend.size();j++)
    {
    for (k=posbeg[j];k<=posend[j];k++,workindex++)
      ind(*workindex,0) = j;
    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\ind.res");
  // ind.prettyPrint(out);
  // TEST


  // 5. Compute meaneffectnr, mclosest, effectvalues

  double dm_mean1 = dm.mean(0);
  double dm_mean2 = dm.mean(1);
  effectvalues.erase(effectvalues.begin(),effectvalues.end());
  double d1,d2;
  meaneffectnr = 0;
  double distclosest,distcurrent;
  distclosest = pow(data(posbeg[0],0)-dm_mean1,2)+
                pow(data(posbeg[0],1)-dm_mean2,2);

  for(j=0;j<posbeg.size();j++)
    {
    d1 = data(posbeg[j],0);
    d2 = data(posbeg[j],1);
    distcurrent = pow(d1-dm_mean1,2)+pow(d2-dm_mean2,2);
    if ( distcurrent < distclosest)
      {
      meaneffectnr = j;
      distclosest = distcurrent;
      }

    effectvalues.push_back(ST::doubletostring(d1,0) + "  "
                           + ST::doubletostring(d2,0));

    }


  compute_meaneffectintvar();

  }


double DESIGN_kriging::compute_matern(double & nu,double & r)
  {

  if(nu==0.5)
    {
    return exp(-r);
    }
  else if(nu==1.5)
    {
    return exp(-r)*(1+r);
    }
  else if(nu==2.5)
    {
    return exp(-r)*(1+r+r*r/3);
    }
  else
    {
    return exp(-r)*(1+r+2*r*r/5+r*r*r/15);
    }
  }


void DESIGN_kriging::compute_tildeZ(void)
  {
  Zout = datamatrix(posbeg.size(),nrpar,0);

  unsigned i,j;
  double r;

  for(i=0; i<Zout.rows(); i++)
    {
    for(j=0; j<Zout.cols(); j++)
      {
      r=sqrt( pow(xvalues[i]-xknots[j],2) + pow(yvalues[i]-yknots[j],2))/rho;

      Zout(i,j) = compute_matern(nu,r);
      }
    }

  tildeZ_t = Zout.transposed();

  }


void DESIGN_kriging::compute_penalty(void)
  {
  unsigned i,j;
  double r;
  Kfull = datamatrix(xknots.size(),xknots.size(),0);
  for(i=0; i<Kfull.rows(); i++)
    {
    for(j=0; j<Kfull.cols(); j++)
      {
      r=sqrt((xknots[i]-xknots[j])*(xknots[i]-xknots[j]) +
      (yknots[i]-yknots[j])*(yknots[i]-yknots[j]))/rho;

      Kfull(i,j) = compute_matern(nu,r);
      }
    }

  rankK = nrpar;

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\K.res");
  // Kfull.prettyPrint(out);
  // TEST
  }


double DESIGN_kriging::penalty_compute_quadform(datamatrix & beta)
  {
  return Kfull.compute_quadform(beta,0);
  }


void DESIGN_kriging::compute_XtransposedWres(datamatrix & partres, double l)
  {
  XWres.mult(tildeZ_t,partres);
  XWres_p = &XWres;
  }


void DESIGN_kriging::compute_XtransposedWX(void)
  {

  // ofstream out("c:\\bayesx\\testh\\results\\wsum.res");
  // Wsum.prettyPrint(out);

  WsumtildeZ.multdiagfront(Zout,Wsum);
  XWXfull.mult(tildeZ_t,WsumtildeZ);

  // TEST
  ofstream out("c:\\bayesx\\testh\\results\\XWXfull.res");
  XWXfull.prettyPrint(out);

  }


void DESIGN_kriging::compute_precision(double l)
  {

  if (precisiondeclared==false)
    {
    precision = envmatdouble(XWXfull);
    precisiondeclared = true;
    }

  precision.addto(XWX,K,1.0,l);

  }


void DESIGN_kriging::outoptions(GENERAL_OPTIONS * op)
  {

  op->out("  Correlation function: Matern \n");
  op->out("  Parameter nu: " + ST::doubletostring(nu) + "\n");
  op->out("  Parameter rho: " + ST::doubletostring(rho) + "\n");
  op->out("  Number of knots: " + ST::inttostring(nrknots) + "\n");


  }

} // end: namespace MCMC






