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





#include"dagobject.h"

void estimaterun(dagobject & d)
  {

  bool error = false;

  MCMCoptions mopt(
  #if defined(JAVA_OUTPUT_WINDOW)
  d.adminb_p,
  #endif
  d.iterations.getvalue(),d.burnin.getvalue(),d.step.getvalue(),d.logout);
  mopt.set_nrout(d.iterationsprint.getvalue());






  //--------------------- reading dataset information --------------------------

  dataobject * datap = NULL;               // pointer to dataset

  int objpos = findstatobject(*(d.statobj),d.udata.getusingtext(),"dataset");
  statobject * s;
  if (objpos >= 0)
    {
    s = d.statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      d.outerror("ERROR: " + d.udata.getusingtext() + " is not existing\n");
    else
      d.outerror("ERROR: " + d.udata.getusingtext() + " is not a dataset object\n");
    error = true;
    }
    //------------------ end: reading dataset information ------------------------








  //---------------- reading data, creating designmatrices ---------------------

  datamatrix  D;
  vector<ST::string> modelvarnamesv;

  if (error == false)
  {
	  modelvarnamesv = d.mod.getModelVarnamesAsVector();

	  ST::string ifexpression = d.methods[0].getexpression();

	  datap->makematrix(modelvarnamesv,D,ifexpression);
	  d.errormessages = datap->geterrormessages();

	  if (!d.errormessages.empty())
		  error= true;
    } // end: if error == false

  //------------------------- end: reading data --------------------------------






	unsigned i,j;
	unsigned cc,dd;
	unsigned number_of_discrete = 0;
	unsigned number_of_continuous = 0;

	ST::string title;
	ST::string path;

    bool mixed_case=false;


	// ****** different types of possible objects
	vector<FULLCOND_dag> dags;
	vector<FULLCOND_dag_d> dags_d;
	vector<FULLCOND_dag_ia> dags_ia;
	vector<FULLCOND_dag_ia_mixed> dags_ia_m;

	datamatrix ad(D.cols(),1,0);
	MCMC::IA interact(D);











	//******* different options that can be chosen ***********************************/

	ST::string prior_sig = d.prior_sig.getvalue();			// dag: variance of prior
	ST::string res_file = d.res_file.getvalue();			// rj: path for aggregated results
	ST::string fix_file = d.fix_file.getvalue();			// rj: path with the conditions on d
	ST::string family = d.family.getvalue();				// rj: kind of variable and wether there are interactions or not
	ST::string print_models = d.print_models.getvalue();	// rj: criterion for number of models in the output
	ST::string switch_typ = d.switch_typ.getvalue();		// rj: kind of switch type


	unsigned ty = d.typ.getvalue();					// rj: type of the adjacency matrix to start with
	unsigned int number = d.number.getvalue();		// rj: number of dags/models in the output

	bool print_dags = d.print_dags.getvalue();		// dag: detailed output about the regression models
	bool detail_ia = d.detail_ia.getvalue();		    // dag: detailed output about the interactions

	double value_a = d.value_a.getvalue();			// dag: a of IG(a,b)
	double value_b = d.value_b.getvalue();			// dag: b of IG(a,b)
	double alpha = d.alpha.getvalue();				// rj: criterion for number of models in the output

	vector<char> var_type;							//




















  if (error==false)
  {

    //***** test if the differnt columns if D are binary or continuous ************

	char type;
    for(i=0;i<D.cols();i++)
    {
		j=0;
		type = 'd';
		while(j<D.rows()-1)
		{
			if(D(j,i)!=0 && D(j,i)!=1)
			{
				type = 'c';
				j=D.rows();
				number_of_continuous++;
			}
			j++;
       }
       var_type.push_back(type);
    }

    number_of_discrete = D.cols()-number_of_continuous;

    if(number_of_discrete>0 &&  number_of_continuous>0)
		mixed_case=true;

    if(family != "unspecified")
    {
       if(		(mixed_case==true && family != "mixed")
           ||  (mixed_case==false && family=="mixed")
           ||	(( family=="discrete_ia"||family=="discrete")
				   && number_of_discrete != D.cols())
           ||  (family=="continuous" &&number_of_continuous != D.cols())
       )
	   {
           d.outerror("ERROR: Chosen type of family does not correspond to the data\n");
		   error=true;
	   }
    }

  }//error==false



  if(error==false)
  {
    //specifying the family
    if( family=="unspecified")
    {
        if(mixed_case==true)
            family="mixed";
        else if(number_of_discrete>0)
            family="discrete";
        else
            family= "continuous";
    }


    if(family=="mixed")
        dags_ia_m.reserve( D.cols());
    else if (family == "discrete_ia")
        dags_ia.reserve(D.cols());
    else if (family=="discrete")
        dags_d.reserve(number_of_discrete);
    else
        dags.reserve(number_of_continuous);


    if(detail_ia==true)
    {
      if(family!="discrete_ia")
      {
          d.out("NOTE: The option 'detail_ia' presupposes the option family='discrete_ia'.\n");
      }
    }













    //******** create vectors with FULLCOND_dag_* Objects *************************

    for(i=0,cc=0,dd=0;i<D.cols();i++)
    {
       title = "dag_" + modelvarnamesv[i];
#if defined(__BUILDING_LINUX)
       path = d.defaultpath.to_bstr() + "/temp/" + d.name.to_bstr() + ST::inttostring(i) + ".raw";
#else
       path = d.defaultpath.to_bstr() + "\\temp\\" + d.name.to_bstr() + ST::inttostring(i) + ".raw";
#endif
       if (family == "continuous")       // continuous
       {
            dags.push_back(FULLCOND_dag( value_a, value_b, prior_sig, print_dags, D.getCol(i),
                                   1, i, &mopt, D,title, D.cols(), 1, path));
			dags[cc].setflags(MCMC::norelchange);
            cc++;
       }
       else if (family == "discrete")      // binary without interaction
       {
            dags_d.push_back(FULLCOND_dag_d( value_a, value_b, prior_sig, print_dags, D.getCol(i),
                                   1, i, &mopt, D,title, D.cols(), 1, path));
			dags_d[dd].setflags(MCMC::norelchange);
            dd++;
       }
       else if (family == "discrete_ia")     // binary with interactions or mixed with interactions
       {
            dags_ia.push_back(FULLCOND_dag_ia( detail_ia, &interact, value_a, value_b, prior_sig, print_dags, D.getCol(i),
                                    1, i, &mopt, D,title, D.cols(), 1, path));
            dags_ia[i].setflags(MCMC::norelchange);
       }
       else if(family == "mixed")
       {
            dags_ia_m.push_back(FULLCOND_dag_ia_mixed( detail_ia, &interact, value_a, value_b, prior_sig, print_dags,
                                                       D.getCol(i), 1, i, &mopt, D,title, D.cols(), 1, path));
            dags_ia_m[i].setflags(MCMC::norelchange);
       }
       else
           d.outerror("ERROR!!!\n");

    }





    //***********	create vectors with pointers to FULLCOND_dag_* objects		***********
    //***********	one vector contains only pointers to one type of object		***********

    vector<FULLCOND_dag*> dagp;
    vector<FULLCOND_dag_ia*> dagp_ia;
    vector<FULLCOND_dag_ia_mixed*> dagp_ia_m;

	if(family=="discrete_ia")            // with interactions (discrete or mixed case)
    {
		for (i=0;i<D.cols();i++)
			dagp_ia.push_back(&dags_ia[i]);
    }
    else if(mixed_case==false)
    {
        for (i=0,cc=0,dd=0;i<D.cols();i++)
        {
             if(var_type[i]=='d')         // discrete without interacions
             {
                 dagp.push_back(&dags_d[dd]);
                 dd++;
             }
             else if(var_type[i]=='c')    // continuous without interacions
             {
                 dagp.push_back(&dags[cc]);
                 cc++;
             }
        }
   }
   else //mixed_case=0true
   {
        for (i=0;i<D.cols();i++)
			dagp_ia_m.push_back(& (dags_ia_m[i]));
   }







   // datamatrix res(D.rows(),1,1);			// wozu? altes relikt?

#if defined(__BUILDING_LINUX)
   path = d.defaultpath.to_bstr() + "/temp/" + d.name.to_bstr() + "rj.raw";
#else
   path = d.defaultpath.to_bstr() + "\\temp\\" + d.name.to_bstr() + "rj.raw";
#endif



	//***********	create suitable FULLCOND_rj_* object		***********
    //***********	and put it togehter with dags in vector fc	***********

    FULLCOND_rj RJ;
    FULLCOND_rj_int RJ_ia;
    FULLCOND_rj_mix RJ_mix;

	vector<FULLCOND *> fc;

    if (family != "discrete_ia" && mixed_case==false)        // no interactions, discrete or continuous
    {
         RJ = FULLCOND_rj (fix_file, res_file, number, alpha, switch_typ, print_models, ty,dagp,
                                            &mopt,D,"dag_rj",D.cols(),D.cols(),path);
         RJ.setflags(MCMC::norelchange);
         fc.push_back(&RJ);
         for (cc=0, dd=0,i=0;i<D.cols();i++)
         {
             if(var_type[i]=='c')
             {
                 fc.push_back(&dags[cc]);
                 cc++;
             }
             else
             {
                 fc.push_back(&dags_d[dd]);
                 dd++;
             }
         }
    }
    else if   (family == "discrete_ia")  //interactions, discrete
    {
         RJ_ia = FULLCOND_rj_int (fix_file, res_file, number, alpha, switch_typ, print_models, ty, dagp_ia,
                                                 &mopt,D,"dag_rj",D.cols(),D.cols(),path);
         RJ_ia.setflags(MCMC::norelchange);
         fc.push_back(&RJ_ia);
         for (i=0;i<D.cols();i++)
             fc.push_back(&dags_ia[i]);
    }
    else if (mixed_case==true)
    {
         RJ_mix = FULLCOND_rj_mix (fix_file, res_file, number, alpha, switch_typ, print_models, ty, dagp_ia_m,
                                                 &mopt,D,"dag_rj",D.cols(),D.cols(),path);
         RJ_mix.setflags(MCMC::norelchange);
         fc.push_back(&RJ_mix);
         for (i=0;i<D.cols();i++)
             fc.push_back(&dags_ia_m[i]);
	}





	//******* let simulation run ***************************************************

    MCMCsimulate sim(&mopt,fc);
    vector<ST::string> header;
    header.push_back("DAG object: method estimate");
    sim.simulate(header,-1,false);

    } // end: if (error==false)

  }










void dagobject::create(void)
{

	// for method estimate

	// SYNTAX OF COMMANDS:
	// name [model] [weight varname] [by varname] [if expression]
	//      [, options] [using usingtext]


	mod = modelStandard();
	udata = use();

	iterations		= intoption("iterations",52000,1,10000000);
	burnin			= intoption("burnin",2000,0,500000);
	step			= intoption("step",50,1,1000);
	iterationsprint = intoption("printit",100,1,100000000);
	typ				= intoption("type",0,0,4);
	number			= intoption("number",10,0,10000);

	alpha	= doubleoption("alpha",0.05,0,1);
	value_a = doubleoption("delta",1,0,20);
	value_b = doubleoption("lambda",0.005,0,20);

	print_dags = simpleoption("print_dags",false);
	detail_ia = simpleoption("detail_ia",false);

	res_file	= stroption("file_of_results");
	fix_file	= stroption("fix");

	vector<ST::string> priorSig;
	priorSig.push_back("inf");
	priorSig.push_back("non_inf");
	prior_sig = stroption("priori_sigma",priorSig,"non_inf");

	vector<ST::string> kind_var;
    kind_var.push_back("unspecified");
	kind_var.push_back("continuous");
	kind_var.push_back("discrete");
	kind_var.push_back("discrete_ia");
	kind_var.push_back("mixed");
	family = stroption("family",kind_var,"unspecified");

	vector<ST::string> print_mod;
	print_mod.push_back("all");
	print_mod.push_back("prob");
	print_mod.push_back("limit");
	print_mod.push_back("normal");
	print_models = stroption("print_models",print_mod,"normal");

	vector<ST::string> switch_t;
	switch_t.push_back("normal");
	switch_t.push_back("equi");
	switch_t.push_back("mix");
	switch_typ = stroption("switch",switch_t,"normal");


	estimateoptions.push_back(&iterations);
	estimateoptions.push_back(&burnin);
	estimateoptions.push_back(&step);
	estimateoptions.push_back(&iterationsprint);
	estimateoptions.push_back(&typ);
	estimateoptions.push_back(&print_models);
	estimateoptions.push_back(&switch_typ);
	estimateoptions.push_back(&number);
	estimateoptions.push_back(&alpha);
	estimateoptions.push_back(&print_dags);
	estimateoptions.push_back(&res_file);
	estimateoptions.push_back(&fix_file);
	estimateoptions.push_back(&prior_sig);
	estimateoptions.push_back(&value_a);
	estimateoptions.push_back(&value_b);
	estimateoptions.push_back(&family);
	estimateoptions.push_back(&detail_ia);

	// SYNTAX OF COMMANDS:
	// name [model] [weight varname] [by varname] [if expression]
	//      [, options] [using usingtext]

	methods.push_back(command("estimate",&mod,&estimateoptions,&udata,required,
					notallowed,notallowed,optional,optional,required));

	functions[0] = estimaterun;

  }












dagobject::dagobject(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  const ST::string & n,ofstream * lo,istream * i,ST::string p,
            vector<statobject*> * st)
         : statobject(
         #if defined(JAVA_OUTPUT_WINDOW)
         adb,
         #endif
         n,"dag",lo,i,p)
  {

  statobj = st;
  create();

  }


dagobject::dagobject(const dagobject & d) : statobject(statobject(d))
  {
  create();
  statobj = d.statobj;
  }



const dagobject & dagobject::operator=(const dagobject & d)
  {
  create();
  statobj = d.statobj;
  if (this == & d)
	 return *this;
  statobject::operator=(statobject(d));
  return *this;
  }



int dagobject::parse(const ST::string & c)
  {
  optionlist globaloptions = optionlist();
  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

 return(pos);
 }




