// TO FIX: any way to merge this file with
// all_kernels.cl?


// TO FIX: is this needed???
// con2prim_MHD_ should only be called
// if OpenCL is turned off
#if OPENCL
  #include "nuc_eos_cl.h"
#else
  #include "nuc_eos_c.h"
#endif

//# include "NR.cl"
//# include "nuc_eos.cl"
//# include "findtemp.c"
//# include "linterp_for_temp.cl"
//# include "linterp_some.cl"
//# include "readtable.c"

extern void calc_prim(int EOS, double * EOS_const, double S_squared, double BdotS, double B_squared, double * con, double * prim, 
		      double * x, double B1_cov, double B2_cov, double B3_cov,struct nuc_eos_vars *struct_ptr);

extern int bisectionmethod(int EOS, double * EOS_const, double S_squared, double BdotS, double B_squared, double * con, double * prim, 
		    double * x, double B1_cov, double B2_cov, double B3_cov, int CODETEST, struct nuc_eos_vars *struct_ptr);

extern int NR_3D(int EOS, double * EOS_const, double S_squared, double BdotS, double B_squared, double * con, double * prim, 
	  double * x, double g_con[4][4], double B1_cov, double B2_cov, double B3_cov, int CODETEST, struct nuc_eos_vars *struct_ptr, int stepsize);

int con2prim_MHD_(int EOS, double * EOS_const, 
                  double * prim, double * con, 
                  double g_con[4][4], double g_cov[4][4], double g_det, 
                  int CODETEST,
                  int nrho_in, int ntemp_in, int nye_in,
                  double *fourtables_in,     double *logrho_in, 
                  double *logtemp_in,        double *yes_in,
                  double energy_shift_in,    double dtemp_in, 
                  double dtempi_in,          double drho_in, 
                  double drhoi_in,           double dye_in, 
                  double dyei_in, 
                  double eos_rhomax_in,      double eos_rhomin_in, 
                  double eos_tempmin_in,     double eos_tempmax_in, 
                  double eos_yemin_in,       double eos_yemax_in, 
                  int *ivs_short_in){

  struct nuc_eos_vars myvars;
  if(EOS==3){
    myvars.nrho = nrho_in;
    myvars.ntemp = ntemp_in;
    myvars.nye = nye_in;
    myvars.fourtables = fourtables_in;
    
    // debugging
    //printf("\n\n\n\n\nTABLE: %e %e %e\n",myvars.fourtables[0],myvars.fourtables[myvars.nrho],myvars.fourtables[myvars.nrho*myvars.ntemp]);
    
    myvars.logrho = logrho_in;
    myvars.logtemp = logtemp_in;
    myvars.yes = yes_in;
    myvars.energy_shift = energy_shift_in;
    myvars.dtemp = dtemp_in;
    myvars.dtempi = dtempi_in;
    myvars.drho = drho_in;
    myvars.drhoi = drhoi_in;
    myvars.dye = dye_in;
    myvars.dyei = dyei_in;
    myvars.eos_rhomax = eos_rhomax_in;
    myvars.eos_rhomin = eos_rhomin_in;
    myvars.eos_tempmin = eos_tempmin_in;
    myvars.eos_tempmax = eos_tempmax_in;
    myvars.eos_yemin = eos_yemin_in;
    myvars.eos_yemax = eos_yemax_in;
    myvars.ivs_short = ivs_short_in;
  }
 

  // Lower indices - covariant
  double B1_cov = g_cov[1][1]*con[B1_con]+g_cov[1][2]*con[B2_con]+g_cov[1][3]*con[B3_con];
  double B2_cov = g_cov[1][2]*con[B1_con]+g_cov[2][2]*con[B2_con]+g_cov[2][3]*con[B3_con];
  double B3_cov = g_cov[1][3]*con[B1_con]+g_cov[2][3]*con[B2_con]+g_cov[3][3]*con[B3_con];
 
  // Need to calculate for (84) & (85)
  double BdotS = con[B1_con]*con[S1_cov] + con[B2_con]*con[S2_cov] + con[B3_con]*con[S3_cov];
  double B_squared = con[B1_con]*B1_cov + con[B2_con]*B2_cov + con[B3_con]*B3_cov;
  
  // Raise indices
  double S1_con = g_con[1][1]*con[S1_cov]+g_con[1][2]*con[S2_cov]+g_con[1][3]*con[S3_cov];
  double S2_con = g_con[1][2]*con[S1_cov]+g_con[2][2]*con[S2_cov]+g_con[2][3]*con[S3_cov];
  double S3_con = g_con[1][3]*con[S1_cov]+g_con[2][3]*con[S2_cov]+g_con[3][3]*con[S3_cov];

  // S^2 = S^i * S_i
  double S_squared = con[S1_cov]*S1_con + con[S2_cov]*S2_con + con[S3_cov]*S3_con;
  
  double x[3];
  x[0] = 0.0;
  x[1] = 0.0;
  x[2] = 0.0;
  
  // Calculate W, Z
  int count;
  int stepsize;
  stepsize=1; 
  if(EOS==1)
    count = bisectionmethod(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x,B1_cov, B2_cov,B3_cov,CODETEST,&myvars);
  else
    count = NR_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, g_con, B1_cov, B2_cov, B3_cov,CODETEST,&myvars,stepsize);

  if (count > 127) {
    count=0;
    stepsize=2; 
    printf("Inside C2P, now calling with reduced step-size: count = %d, stepsize= %d \n",count,stepsize);
    count = NR_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, g_con, B1_cov, B2_cov, B3_cov,CODETEST,&myvars,stepsize);
    printf("Inside C2P, after calling with reduced step-size: count = %d, stepsize= %d \n",count,stepsize);
    stepsize=1; 
  }
  if (count > 127) {
    count=0;
    stepsize=3; 
    printf("Inside C2P, now calling with reduced step-size: count = %d, stepsize= %d \n",count,stepsize);
    count = NR_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, g_con, B1_cov, B2_cov, B3_cov,CODETEST,&myvars,stepsize);
    printf("Inside C2P, after calling with reduced step-size: count = %d, stepsize= %d \n",count,stepsize);
    stepsize=1; 
  }
  if (count > 127) {
    count=0;
    stepsize=1;
    CODETEST=1;
    printf("Inside C2P, now calling with safe-guess initial values: CODETEST = %d count = %d, stepsize= %d \n",CODETEST,count,stepsize);
    count = NR_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, g_con, B1_cov, B2_cov, B3_cov,CODETEST,&myvars,stepsize);
    printf("Inside C2P, after calling with safe-guess initial values: CODETEST = %d count = %d, stepsize= %d \n",CODETEST,count,stepsize);
    CODETEST=0;
    stepsize=1;
  }
  if (count > 127) {
    count=0;
    stepsize=2;
    CODETEST=1;
    printf("Inside C2P, now calling with reduced step-size and safe-guess initial values: CODETEST = %d count = %d, stepsize= %d \n",CODETEST,count,stepsize);
    count = NR_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, g_con, B1_cov, B2_cov, B3_cov,CODETEST,&myvars,stepsize);
    printf("Inside C2P, after calling with reduced step-sizeand safe-guess initial values: CODETEST = %d count = %d, stepsize= %d \n",CODETEST,count,stepsize);
    CODETEST=0;
    stepsize=1;
  }
  if (count > 127) {
    count=0;
    stepsize=5;
    CODETEST=0;
    printf("Inside C2P, now calling with reduced step-size and safe-guess initial values: CODETEST = %d count = %d, stepsize= %d \n",CODETEST,count,stepsize);
    count = NR_3D(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, g_con, B1_cov, B2_cov, B3_cov,CODETEST,&myvars,stepsize);
    printf("Inside C2P, after calling with reduced step-sizeand safe-guess initial values: CODETEST = %d count = %d, stepsize= %d \n",CODETEST,count,stepsize);
    CODETEST=0;
    stepsize=1;
  }

 
  calc_prim(EOS, EOS_const, S_squared, BdotS, B_squared, con, prim, x, B1_cov, B2_cov, B3_cov,&myvars);
  return count;
}


