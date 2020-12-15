double linterp2D(__global double* xs, __global double* ys, double* fs, double x, double y);

void bisection(double lr, double lt0, double ye, double leps0, double* ltout,
	       int iv, double prec, int* keyerrt, struct nuc_eos_vars *struct_ptr);

void nuc_eosK_C_findtemp(double lr, double lt0, double ye, 
			double lepsin, double prec, double *ltout,
			int *keyerrt,struct nuc_eos_vars *struct_ptr) {

  // local vars
  int itmax = 10; // use at most 10 iterations, then go to bisection
  double dlepsdlt; // derivative dlogeps/dlogT
  double ldt;
  double leps,leps0,leps1; // temp vars for eps
  double ltn, lt, lt1; // temp vars for temperature
  double ltmax = struct_ptr->logtemp[struct_ptr->ntemp-1]; // max temp
  double ltmin = struct_ptr->logtemp[0]; // min temp

  // setting up some vars
  *keyerrt = 0;
  leps0 = lepsin;
  leps1 = leps0;
  lt = lt0;
  lt1 = lt;

  // step 1: do we already have the right temperature
  nuc_eosK_C_linterp_for_temp(lr,lt,ye,&leps,struct_ptr->fourtables,struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,
			     struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,&dlepsdlt);

  if(fabs(leps-leps0) < prec*fabs(leps0)) {
    *ltout = lt0;
    return;
  }
  lt1 = lt;
  leps1 = leps;

  int it = 0;
  while(it < itmax) {

    // step 2: check if the two bounding values of the temperature
    //         give eps values that enclose the new eps.
    int itemp = MIN(MAX(1 + (int)(( lt - struct_ptr->logtemp[0] - 1.0e-10) * struct_ptr->dtempi),1),struct_ptr->ntemp-1); 
    int irho = MIN(MAX(1 + (int)(( lr - struct_ptr->logrho[0] - 1.0e-10) * struct_ptr->drhoi),1),struct_ptr->nrho-1); 
    int iye = MIN(MAX(1 + (int)(( ye - struct_ptr->yes[0] - 1.0e-10) * struct_ptr->dyei),1),struct_ptr->nye-1); 

    double epst1, epst2;
    // lower temperature
    //{
      // get data at 4 points
      double fs[4];
      // point 0
      int ifs = 1 + NSUBTABLES*(irho-1 + struct_ptr->nrho*((itemp-1) + struct_ptr->ntemp*(iye-1)));
      fs[0] = struct_ptr->fourtables[ifs];
      // point 1 
      ifs = 1 + NSUBTABLES*(irho + struct_ptr->nrho*((itemp-1) + struct_ptr->ntemp*(iye-1)));
      fs[1] = struct_ptr->fourtables[ifs];
      // point 2 
      ifs = 1 + NSUBTABLES*(irho-1 + struct_ptr->nrho*((itemp-1) + struct_ptr->ntemp*(iye)));
      fs[2] = struct_ptr->fourtables[ifs];
      // point 3
      ifs = 1 + NSUBTABLES*(irho + struct_ptr->nrho*((itemp-1) + struct_ptr->ntemp*(iye)));
      fs[3] = struct_ptr->fourtables[ifs];
      
      epst1 = linterp2D(&struct_ptr->logrho[irho-1],&struct_ptr->yes[iye-1], fs, lr, ye);
    //}
    // upper temperature
    //{
      // get data at 4 points
      //double fs[4];
      // point 0
      ifs = 1 + NSUBTABLES*(irho-1 + struct_ptr->nrho*((itemp) + struct_ptr->ntemp*(iye-1)));
      fs[0] = struct_ptr->fourtables[ifs];
      // point 1 
      ifs = 1 + NSUBTABLES*(irho + struct_ptr->nrho*((itemp) + struct_ptr->ntemp*(iye-1)));
      fs[1] = struct_ptr->fourtables[ifs];
      // point 2 
      ifs = 1 + NSUBTABLES*(irho-1 + struct_ptr->nrho*((itemp) + struct_ptr->ntemp*(iye)));
      fs[2] = struct_ptr->fourtables[ifs];
      // point 3
      ifs = 1 + NSUBTABLES*(irho + struct_ptr->nrho*((itemp) + struct_ptr->ntemp*(iye)));
      fs[3] = struct_ptr->fourtables[ifs];
      
      epst2 = linterp2D(&struct_ptr->logrho[irho-1],&struct_ptr->yes[iye-1], fs, lr, ye);
    //}
    
    // Check if we are already bracketing the input internal
    // energy. If so, interpolate for new T.
    if(leps0 >= epst1 && leps0 <= epst2) {
      
      *ltout = (struct_ptr->logtemp[itemp]-struct_ptr->logtemp[itemp-1]) / (epst2 - epst1) * 
	(leps0 - epst1) + struct_ptr->logtemp[itemp-1];
      return;

    }

    // well, then do a Newton-Raphson step
    ldt = -(leps - leps0) / dlepsdlt;
    ltn = MIN(MAX(lt + ldt,ltmin),ltmax);
    lt1 = lt;
    lt = ltn;
    leps1 = leps;

    nuc_eosK_C_linterp_for_temp(lr,lt,ye,&leps,struct_ptr->fourtables,struct_ptr->nrho,struct_ptr->ntemp,struct_ptr->nye,
			       struct_ptr->logrho,struct_ptr->logtemp,struct_ptr->yes,&dlepsdlt);

    if(fabs(leps-leps0) < prec*fabs(leps0)) {
      *ltout = lt;
      return;
    }
    
    // if we are closer than 10^-3  to the 
    // root (eps-eps0)=0, we are switching to 
    // the secant method, since the table is rather coarse and the
    // derivatives may be garbage.
    if(fabs(leps-leps0) < 1.0e-3*fabs(leps0)) {
      dlepsdlt = (leps-leps1)/(lt-lt1);
    }

    it++;
  }

  if(it >= itmax - 1) {
    // try bisection

    bisection(lr, lt0, ye, leps0, ltout, 1, prec, keyerrt,struct_ptr);

    return;
  }



}

double linterp2D(__global double* xs, __global double* ys, double* fs, double x, double y)
{

  //  2     3 
  //
  //  0     1
  //
  // first interpolate in x between 0 and 1, 2 and 3
  // then interpolate in y
  // assume rectangular grid
  
  double t1 = (fs[1]-fs[0])/(xs[1]-xs[0]) * (x - xs[0]) + fs[0];
  double t2 = (fs[3]-fs[2])/(xs[1]-xs[0]) * (x - xs[0]) + fs[2];

  return (t2 - t1)/(ys[1]-ys[0]) * (y-ys[0]) + t1;
}

void bisection(double lr, double lt0, double ye, double leps0, double* ltout,
	       int iv, double prec, int* keyerrt, struct nuc_eos_vars *struct_ptr) {

  // iv is the index of the table variable we do the bisection on

  int bcount = 0; 
  int maxbcount = 80;
  int itmax = 50;
  
  // temporary local vars
  double lt, lt1, lt2;
  double ltmin = struct_ptr->logtemp[0];
  double ltmax = struct_ptr->logtemp[struct_ptr->ntemp-1];
  double f1,f2,fmid,dlt,ltmid;
  double f1a = 0.0;
  double f2a = 0.0;

  // prepare
  lt = lt0;
  lt1 = log10( MIN(pow(10.0,ltmax),1.1*(pow(10.0,lt0))) );
  lt2 = log10( MAX(pow(10.0,ltmin),0.9*(pow(10.0,lt0))) );

  int nvars = 1;
  int ivs[1] = {iv};
  nuc_eosK_C_linterp_one(lr, lt1, ye, &f1a, struct_ptr->fourtables, 
  			 ivs, struct_ptr->nrho, struct_ptr->ntemp, struct_ptr->nye,
  			 struct_ptr->logrho, struct_ptr->logtemp, struct_ptr->yes, struct_ptr);

  nuc_eosK_C_linterp_one(lr, lt2, ye, &f2a, struct_ptr->fourtables, 
			 ivs, struct_ptr->nrho, struct_ptr->ntemp, struct_ptr->nye,
  			 struct_ptr->logrho, struct_ptr->logtemp, struct_ptr->yes, struct_ptr);


  f1=f1a-leps0;
  f2=f2a-leps0;
  
  // iterate until we bracket the right eps, but enforce
  // dE/dt > 0, so eps(lt1) > eps(lt2)
#if 0
  int ifixdeg = 0;
  int ifixdeg_max = 20;
#endif
  while(f1*f2 >= 0.0) {
    
    lt1 = log10( MIN(pow(10.0,ltmax),1.2*(pow(10.0,lt1))) );
    lt2 = log10( MAX(pow(10.0,ltmin),0.8*(pow(10.0,lt2))) );
    nuc_eosK_C_linterp_one(lr, lt1, ye, &f1a, struct_ptr->fourtables, 
			   ivs, struct_ptr->nrho, struct_ptr->ntemp, struct_ptr->nye,
			   struct_ptr->logrho, struct_ptr->logtemp, struct_ptr->yes,struct_ptr);
    nuc_eosK_C_linterp_one(lr, lt2, ye, &f2a, struct_ptr->fourtables, 
			   ivs, struct_ptr->nrho, struct_ptr->ntemp, struct_ptr->nye,
			   struct_ptr->logrho, struct_ptr->logtemp, struct_ptr->yes,struct_ptr);

#if 0
    // special enforcement of eps(lt1)>eps(lt2)
    while(f1a < f2a && ifixdeg < ifixdeg_max) {
      lt1 = log10( MIN(pow(10.0,ltmax),1.2*(pow(10.0,lt1))) );
      nuc_eosK_C_linterp_one(lr, lt1, ye, &f1a, fourtables, 
			     ivs, nrho, ntemp, nye,
			     logrho, logtemp, yes);
      ifixdeg++;

    }
#endif


    f1=f1a-leps0;
    f2=f2a-leps0;





    bcount++;
    if(bcount >= maxbcount) {
      *keyerrt = 667;
      return;
    }

  }

  if(f1 < 0.0) {
    lt = lt1;
    dlt = log10( pow(10.0,lt2) - pow(10.0,lt1) );
  } else {
    lt = lt2;
    dlt = log10( pow(10.0,lt1) - pow(10.0,lt2) );
  }


  int it;
  for(it=0;it<itmax;it++) {
    dlt = log10( pow(10.0,dlt) * 0.5 );
    ltmid = log10( pow(10.0,lt) + pow(10.0,dlt) );
    nuc_eosK_C_linterp_one(lr, ltmid, ye, &f2a, struct_ptr->fourtables, 
			   ivs, struct_ptr->nrho, struct_ptr->ntemp, struct_ptr->nye,
			   struct_ptr->logrho, struct_ptr->logtemp, struct_ptr->yes,struct_ptr);

    fmid=f2a-leps0;
    if(fmid <= 0.0) lt=ltmid;


    if(fabs(1.0-f2a/leps0) <= prec) {
      *ltout = ltmid;
      return;
    }

  }

  if(it >= itmax-1) {
    *keyerrt = 667;
    return;
  }

}
