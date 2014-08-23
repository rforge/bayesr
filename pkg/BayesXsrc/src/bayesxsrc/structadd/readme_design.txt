  // FUNCTION: init_data
  // TASK: sorts the data such that the precision has minimum envelope
  //       computes index_data
  //       computes Zout, posbeg, posend
  //       computes nrpar
  //       computes effectvalues
  //       computes consecutive
  //       computes ZoutT, index_ZoutT    


  virtual void init_data(const datamatrix & dm, const datamatrix & iv);


  // FUNCTION: read_options
  // TASK: reads options and initializes varnames stored in datanames

  virtual void read_options(vector<ST::string> & op,vector<ST::string> & vn);


  // FUNCTION: compute_penalty
  // TASK: computes the penalty matrix and determines rankK

  virtual void compute_penalty(void);



  OPTIONAL: 


  // FUNCTION: computes XWres
  // TASK: computes XWres, res is the partial residual

  virtual void compute_XtransposedWres(datamatrix & partres, double l);

  // FUNCTION: compute_XtransposedWX_XtransposedWres
  // TASK: computes XWX and XWres, res is the partial residual

  virtual void compute_XtransposedWX(void);

  // FUNCTION: compute_XtransposedWX_XtransposedWres
  // TASK: computes XWX and XWres, res is the partial residual

  virtual void compute_XtransposedWX_XtransposedWres(datamatrix & partres, double l);




